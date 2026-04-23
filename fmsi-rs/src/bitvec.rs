use std::io::{Read, Write};

#[derive(Clone, Debug, Default)]
pub struct BitVecRank {
    len: usize,
    words: Vec<u64>,
    checkpoints: Vec<u64>,
}

impl BitVecRank {
    pub fn with_len(len: usize) -> Self {
        let words = vec![0u64; len.div_ceil(64)];
        Self {
            len,
            words,
            checkpoints: Vec::new(),
        }
    }

    pub fn from_bools<I: IntoIterator<Item = bool>>(iter: I) -> Self {
        let vals: Vec<bool> = iter.into_iter().collect();
        let mut bv = Self::with_len(vals.len());
        for (i, bit) in vals.iter().copied().enumerate() {
            bv.set(i, bit);
        }
        bv.rebuild_rank();
        bv
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn set(&mut self, idx: usize, value: bool) {
        assert!(idx < self.len);
        let word = idx / 64;
        let bit = idx % 64;
        if value {
            self.words[word] |= 1u64 << bit;
        } else {
            self.words[word] &= !(1u64 << bit);
        }
    }

    pub fn get(&self, idx: usize) -> bool {
        assert!(idx < self.len);
        let word = idx / 64;
        let bit = idx % 64;
        ((self.words[word] >> bit) & 1) != 0
    }

    pub fn rebuild_rank(&mut self) {
        self.checkpoints = Vec::with_capacity(self.words.len() + 1);
        let mut prefix = 0u64;
        self.checkpoints.push(prefix);
        for &word in &self.words {
            prefix += word.count_ones() as u64;
            self.checkpoints.push(prefix);
        }
    }

    pub fn rank1(&self, pos: usize) -> usize {
        assert!(pos <= self.len);
        let word = pos / 64;
        let bit = pos % 64;
        let mut total = self.checkpoints[word] as usize;
        if bit > 0 && word < self.words.len() {
            let mask = (1u64 << bit) - 1;
            total += (self.words[word] & mask).count_ones() as usize;
        }
        total
    }

    pub fn rank0(&self, pos: usize) -> usize {
        pos - self.rank1(pos)
    }

    pub fn write_to<W: Write>(&self, mut w: W) -> std::io::Result<()> {
        w.write_all(b"BV64")?;
        w.write_all(&(self.len as u64).to_le_bytes())?;
        w.write_all(&(self.words.len() as u64).to_le_bytes())?;
        for word in &self.words {
            w.write_all(&word.to_le_bytes())?;
        }
        Ok(())
    }

    pub fn read_from<R: Read>(mut r: R) -> std::io::Result<Self> {
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if &magic != b"BV64" {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "invalid bitvector magic",
            ));
        }
        let mut buf = [0u8; 8];
        r.read_exact(&mut buf)?;
        let len = u64::from_le_bytes(buf) as usize;
        r.read_exact(&mut buf)?;
        let words_len = u64::from_le_bytes(buf) as usize;
        let mut words = vec![0u64; words_len];
        for word in &mut words {
            r.read_exact(&mut buf)?;
            *word = u64::from_le_bytes(buf);
        }
        let mut out = Self {
            len,
            words,
            checkpoints: Vec::new(),
        };
        out.rebuild_rank();
        Ok(out)
    }
}
