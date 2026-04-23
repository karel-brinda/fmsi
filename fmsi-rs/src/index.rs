use crate::bitvec::BitVecRank;
use crate::dna::is_upper;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug, Default)]
pub struct StrandPredictor {
    pub score: i32,
    pub result_scores: [i32; 2],
    pub previous: Option<usize>,
}

impl StrandPredictor {
    #[inline]
    fn clipped(x: i32, clipper: i32) -> i32 {
        x.max(-clipper).min(clipper)
    }

    pub fn log_result(&mut self, forward_query_result: i32, reverse_query_result: i32) {
        let difference = forward_query_result - reverse_query_result;
        self.score = Self::clipped(self.score + difference, 7);
        if let Some(previous) = self.previous {
            self.result_scores[previous] = Self::clipped(self.result_scores[previous] + difference, 7);
        }
        self.previous = Some((forward_query_result > reverse_query_result) as usize);
    }

    pub fn predict_swap(&self) -> bool {
        if let Some(previous) = self.previous {
            if self.result_scores[previous] != 0 {
                return self.result_scores[previous] < 0;
            }
        }
        self.score < 0
    }
}

#[derive(Clone, Debug, Default)]
pub struct FmsIndex {
    pub ac_gt: BitVecRank,
    pub ac: BitVecRank,
    pub gt: BitVecRank,
    pub sa_transformed_mask: BitVecRank,
    pub counts: [usize; 4],
    pub dollar_position: usize,
    pub klcp: Option<BitVecRank>,
    pub k: usize,
    pub predictor: StrandPredictor,
}

impl FmsIndex {
    pub fn dump_to_prefix<P: AsRef<Path>>(&self, prefix: P) -> Result<(), String> {
        let prefix = prefix.as_ref().to_string_lossy();
        let base = format!("{prefix}.fmsi");

        self.write_bitvec(&self.ac_gt, &(base.clone() + ".ac_gt"))?;
        self.write_bitvec(&self.ac, &(base.clone() + ".ac"))?;
        self.write_bitvec(&self.gt, &(base.clone() + ".gt"))?;
        self.write_bitvec(&self.sa_transformed_mask, &(base.clone() + ".mask"))?;
        if let Some(klcp) = &self.klcp {
            self.write_bitvec(klcp, &(base.clone() + ".klcp"))?;
        }

        let mut w = BufWriter::new(File::create(base + ".misc").map_err(|e| e.to_string())?);
        writeln!(w, "{}", self.dollar_position).map_err(|e| e.to_string())?;
        for c in self.counts {
            writeln!(w, "{}", c).map_err(|e| e.to_string())?;
        }
        writeln!(w, "{}", self.k).map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn load_from_prefix<P: AsRef<Path>>(prefix: P, use_klcp: bool) -> Result<Self, String> {
        let prefix = prefix.as_ref().to_string_lossy();
        let base = format!("{prefix}.fmsi");

        let ac_gt = Self::read_bitvec(&(base.clone() + ".ac_gt"))?;
        let ac = Self::read_bitvec(&(base.clone() + ".ac"))?;
        let gt = Self::read_bitvec(&(base.clone() + ".gt"))?;
        let sa_transformed_mask = Self::read_bitvec(&(base.clone() + ".mask"))?;
        let klcp_path = base.clone() + ".klcp";
        let klcp = if use_klcp && Path::new(&klcp_path).exists() {
            Some(Self::read_bitvec(&klcp_path)?)
        } else {
            None
        };

        let mut reader = BufReader::new(File::open(base + ".misc").map_err(|e| e.to_string())?);
        let mut lines = Vec::new();
        let mut line = String::new();
        while reader.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
            lines.push(line.trim().to_string());
            line.clear();
        }
        if lines.len() < 6 {
            return Err("invalid .misc file".to_string());
        }

        let dollar_position = lines[0].parse::<usize>().map_err(|e| e.to_string())?;
        let counts = [
            lines[1].parse::<usize>().map_err(|e| e.to_string())?,
            lines[2].parse::<usize>().map_err(|e| e.to_string())?,
            lines[3].parse::<usize>().map_err(|e| e.to_string())?,
            lines[4].parse::<usize>().map_err(|e| e.to_string())?,
        ];
        let k = lines[5].parse::<usize>().map_err(|e| e.to_string())?;

        Ok(Self {
            ac_gt,
            ac,
            gt,
            sa_transformed_mask,
            counts,
            dollar_position,
            klcp,
            k,
            predictor: StrandPredictor::default(),
        })
    }

    fn write_bitvec<P: AsRef<Path>>(&self, bv: &BitVecRank, path: P) -> Result<(), String> {
        let mut w = BufWriter::new(File::create(path).map_err(|e| e.to_string())?);
        bv.write_to(&mut w).map_err(|e| e.to_string())
    }

    fn read_bitvec<P: AsRef<Path>>(path: P) -> Result<BitVecRank, String> {
        let mut r = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
        BitVecRank::read_from(&mut r).map_err(|e| e.to_string())
    }

    pub fn export_ms(&self) -> String {
        let masked_letters = b"acgtACGT";
        let n = self.sa_transformed_mask.len() - 1;
        let mut ret = vec![b'N'; n];
        let mut bw_index = 0usize;
        for i in 0..n {
            let letter = crate::search::access(self, bw_index) as usize;
            bw_index = self.counts[letter] + crate::search::rank(self, bw_index, letter as u8);
            let mask_offset = if self.sa_transformed_mask.get(bw_index) { 4 } else { 0 };
            ret[n - 1 - i] = masked_letters[letter + mask_offset];
        }
        String::from_utf8(ret).expect("masked superstring must stay ASCII")
    }

    pub fn from_exported_ms_roundtrip(ms: &str) -> usize {
        ms.as_bytes().iter().filter(|&&c| is_upper(c)).count()
    }
}
