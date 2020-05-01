extern crate rust_htslib;
use std::process;
use std::str::from_utf8;
use std::convert::TryInto;
use std::collections::HashMap;
use std::env;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;


fn main() {
   
    let mut args: Vec<String> = env::args().collect();
	let config = Config::new(&args);

	let r1_bases = process_bam(&config.bam, config.lr1, &config.loci1, config.doqc);
	let r2_bases = process_bam(&config.bam, config.lr2, &config.loci2, config.doqc);
	println!("Genotypes on fragments spanning {}:{} and {}:{} (Insert length {}bp)", &config.chr, &config.loci1, &config.chr, &config.loci2, &config.loci2 - &config.loci1);
	summarize(&r1_bases, &r2_bases);
}

fn summarize_loci(hm: &HashMap<String, String>){
							 	 //A, T, G, C, INDEL
	let mut depth: Vec<f32> = vec![0.0, 0.0, 0.0, 0.0, 0.0];

	for v in hm.values(){
		match v.as_str() {
			"A" => depth[0] = depth[0]+1.0,
			"T" => depth[1] = depth[1]+1.0,
			"G" => depth[2] = depth[2]+1.0,
			"C" => depth[3] = depth[3]+1.0,
			_ => depth[4] = depth[4]+1.0,
		}
	}

	let tdepth = &depth[0] + &depth[1] + &depth[2] + &depth[3] + &depth[4];
	println!("{:?}", tdepth);

	println!("A|T|G|C|INDELS");
	println!("{}({:.2})|{}({:.2})|{}({:.2})|{}({:.2})|{}({:.2})", &depth[0], (&depth[0]/tdepth), 
		&depth[1], &depth[1]/tdepth, &depth[2], &depth[2]/tdepth, &depth[3], &depth[3]/tdepth, &depth[4], &depth[4]/tdepth);

}


fn summarize(hm1: &HashMap<String, String>, hm2: &HashMap<String, String>){

	let mut genotypes: HashMap<String, f32> = HashMap::new();

	for k in hm1.keys(){
		if hm2.contains_key(k){
			let gt: String = hm1[k].clone() + "/" + &hm2[k];
			if genotypes.contains_key(&gt){
				if let Some(x) = genotypes.get_mut(&gt) {
					    *x = *x+1.0;
					}

			}else{
				genotypes.entry(gt).or_insert(1.0);
			}
		}
	}

	let mut tdepth: f32 = 0.0;
	for v in genotypes.values(){
		tdepth = tdepth + v;
	}

	for (k, v) in genotypes{
		println!("{}\t{}({:.2})", k, v, v/tdepth);	
	}
}

fn qc_pass(rec: &Record, mapq_lim: &u8) -> bool{
    if rec.mapq() > *mapq_lim && !rec.is_mate_unmapped() && !rec.is_secondary() && !rec.is_supplementary(){
        true
    }else{
        false
    }
}


fn process_bam(bam_file: &String, region: String, loci: &i64, qc: bool) -> HashMap::<String, String>{
	let mut bam = bam::IndexedReader::from_path(bam_file).unwrap();
	bam.fetch_str(region.as_bytes()).unwrap();
	let mut rinfo: HashMap<String, String> = HashMap::new();

	//					A, T, G, C, Other
	//let mut basecounts: vec![u32] = [0, 0, 0, 0, 0];
	//let mut basecounts: Vec<i32> = vec![0, 0, 0, 0, 0];

    for r in bam.records(){
		let r = r.unwrap();
		//println!("{:?}", r.qname());
		//let base = seq_at(&r, &config.loci1);
		//println!("{}", base);

		if qc == true{
				if qc_pass(&r, &30) == true{
    			let base = seq_at(&r, loci); 
    			rinfo.insert(base.rid, base.base);
    		}
		}else{
			let base = seq_at(&r, loci); 
			rinfo.insert(base.rid, base.base);
    	}
    	
    }

	println!("Base counts at {}", &region);
	summarize_loci(&rinfo);
	println!("----------------------------");
	rinfo
}

fn seq_at(rec: &Record, pos: &i64) -> Rbase {

	let mut pos_on_read: usize = 0;
	match rec.cigar().read_pos((*pos).try_into().unwrap(), false,  true).unwrap(){
			Some(x) => pos_on_read = x as usize,
			None => pos_on_read = 0,
	}

	let seq: String = from_utf8(&rec.seq().as_bytes()).unwrap().to_string();
	let rname: String = from_utf8(&rec.qname()).unwrap().to_string();
	//let mut base_at: String = (&seq[pos_on_read..(pos_on_read+1)]).to_string();
	let mut base_at_manual: String = String::from("NA");
	
	
	let mut temp_pos: i64 = rec.pos().clone();

	for c in &rec.cigar(){
		let cig_len: i64 = c.len().try_into().unwrap();
		

		if c.char() == 'M'{
			temp_pos = temp_pos + cig_len;
			if temp_pos > *pos{
				base_at_manual =  (&seq[pos_on_read..(pos_on_read+1)]).to_string();
				break;
			}
		}

		if c.char() == 'I'{
			if temp_pos == *pos{
				let cig_len: usize = c.len().try_into().unwrap();
				base_at_manual =  (&seq[pos_on_read..(pos_on_read+cig_len)]).to_string();
				break;
			}
		}
	}
	Rbase::new(rname, base_at_manual)
}


#[derive(Debug)]
struct Config {
	bam: String,
	chr: String,
	loci1: i64,
	lr1: String,
	loci2: i64,
	lr2: String,
	doqc: bool,
}

impl Config {
	fn new(args: &[String]) -> Config{

		if args.len() < 4{
			Config::print_usage();
		}

		let qc = "false".to_string();
		if args.len() == 5{
			let qc = args[5].clone();
		}

		let bam = args[1].clone();
		let l1 = args[2].clone();
		let l2 = args[3].clone();
		

		let mut doqc: bool = true;
		match qc.as_str() {
			"true" => doqc = true,
			"false" => doqc = false,
			_ => panic!("Last argument should be true/false"),
		}
		
				
		let l1_spl: Vec<&str> = l1.split(':').collect();
		let l2_spl: Vec<&str> = l2.split(':').collect();
		
		let chr: String = l1_spl[0].parse().unwrap();
		let loci1: i64 =  l1_spl[1].parse().unwrap();
		let loci2: i64 =  l2_spl[1].parse().unwrap();

		let lr1: String = chr.clone() + ":"  + l1_spl[1] + "-" + l1_spl[1];
		let lr2: String = chr.clone() + ":"  + l2_spl[1] + "-" + l2_spl[1];

		if l1_spl[0] != l2_spl[0]{
			println!("Chromosome names must match for both the loci!");
			process::exit(1);
		}else{
			Config{bam, chr, loci1, lr1, loci2, lr2, doqc}	
		}
	}

	fn print_usage(){
		println!("Usage: hapcounter <bam> <region1> <region2> [optionalQC:ture/false]");
		println!("       hapcounter molt4LMO2.bam chr11:33881016 chr11:33956761");
		println!("       hapcounter molt4LMO2.bam chr11:33881016 chr11:33956761 true");
		process::exit(0);
	}
}

#[derive(Debug)]
struct Rbase {
	rid: String,
	base: String,
}

impl Rbase {
	fn new(rid: String, rbase: String) -> Rbase{
		let mut base = String::new();
		match rbase.as_str() {
			"A" | "T" | "G" | "C" => base = rbase,
			_ => base = "INDEL".to_string(),
		}
		Rbase{rid, base}
	}
}