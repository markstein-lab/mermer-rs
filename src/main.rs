// Copyright © 2019-2020 Jakob L. Kreuze <zerodaysfordays@sdf.org>
//
// This file is part of mermer-rs.
//
// mermer-rs is free software; you can redistribute it and/or modify it
// under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// mermer-rs is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with mermer-rs. If not, see <http://www.gnu.org/licenses/>.

#[macro_use]
extern crate anyhow;

use anyhow::Result;
use bio::data_structures::annot_map::AnnotMap;
use bio::io::gff;
use bio::io::gff::Record;
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;
use serde::Serialize;

#[derive(Debug, Serialize)]
struct Sequence {
    #[serde(rename = "sequenceName")]
    name: String,
    #[serde(rename = "resultLocation")]
    location: u64,
    genes: Vec<Gene>,
}

#[derive(Debug, Serialize)]
struct Gene {
    start: u64,
    end: u64,
    #[serde(rename = "geneId")]
    id: String,
    #[serde(rename = "geneName")]
    name: String,
    transcripts: Vec<Transcript>,
}

#[derive(Debug, Serialize)]
struct Transcript {
    start: u64,
    end: u64,
    #[serde(rename = "transcriptId")]
    id: String,
    #[serde(rename = "transcriptName")]
    name: String,
    exons: Vec<Exon>,
    ranges: Vec<Range>,
}

#[derive(Debug, Serialize)]
struct Exon {
    start: u64,
    end: u64,
    #[serde(rename = "exonId")]
    id: String,
    #[serde(rename = "exonNumber")]
    number: u64,
    #[serde(rename = "exonVersion")]
    version: u64,
    ranges: Vec<Range>,
}

#[derive(Debug, Serialize)]
struct Range {
    start: u64,
    end: u64,
    #[serde(rename = "type")]
    range_type: String,
}

fn build_tree<T>(path: T) -> Result<AnnotMap<String, Gene>>
where
    T: AsRef<std::path::Path>,
{
    // This function assumes that the annotation file is strictly ordered such
    // that "gene" features are followed by their "transcript" features, which
    // are then followed by their "exon" features, which are then followed by
    // their "CDS", "UTR", "start_codon", and "stop_codon" features.
    //
    // The naïve solution of using a HashMap for random access of the genes
    // processed so far is slow; I believe that the effort demanded by
    // preprocessing a GTF to fit this assumption is worth the speed benefit.
    let mut annotations = AnnotMap::new();
    let mut reader = gff::Reader::from_file(path, gff::GffType::GFF2)?;
    let (mut gene, mut contig) = (None, None);
    let mut transcript = None;
    let mut exon = None;
    for record in reader.records() {
        let record = record?;
        let attributes = record.attributes();

        match record.feature_type() {
            "gene" => {
                if let (Some(gene), Some(contig)) = (gene, contig) {
                    annotations.insert_at(gene, &contig);
                }
                let id = String::from(attributes.get("gene_id").unwrap());
                let name = String::from(attributes.get("gene_name").unwrap());
                gene = Some(Gene {
                    start: *record.start(),
                    end: *record.end(),
                    id,
                    name,
                    transcripts: Vec::new(),
                });
                contig = Some(Contig::new(
                    record.source().into(),
                    (*record.start()) as isize,
                    (*record.end() - *record.start()) as usize,
                    record.strand().unwrap(),
                ));
            }
            "transcript" => {
                if let Some(transcript) = transcript {
                    gene.as_mut().unwrap().transcripts.push(transcript);
                }

                let id = String::from(attributes.get("transcript_id").unwrap());
                let name = String::from(attributes.get("transcript_name").unwrap());
                transcript = Some(Transcript {
                    start: *record.start(),
                    end: *record.end(),
                    id,
                    name,
                    exons: Vec::new(),
                    ranges: Vec::new(),
                });
            }
            "exon" => {
                if let Some(exon) = exon {
                    transcript.as_mut().unwrap().exons.push(exon);
                }

                let id = String::from(attributes.get("exon_id").unwrap());
                let number = String::from(attributes.get("exon_number").unwrap());
                let version = String::from(attributes.get("exon_version").unwrap());
                exon = Some(Exon {
                    start: *record.start(),
                    end: *record.end(),
                    id,
                    number: number.parse().unwrap(),
                    version: version.parse().unwrap(),
                    ranges: Vec::new(),
                });
            }
            "CDS" | "UTR" | "start_codon" | "stop_codon" => {
                let range = Range {
                    start: *record.start(),
                    end: *record.end(),
                    range_type: String::from(record.feature_type()),
                };

                if attributes.get("exon_id").is_some() {
                    exon.as_mut().unwrap().ranges.push(range);
                } else {
                    transcript.as_mut().unwrap().ranges.push(range);
                }
            }
            _ => {}
        }
    }
    Ok(annotations)
}

fn main() {
    let path = "/home/jakob/University/BIOL 396/dm6/Drosophila_melanogaster.BDGP5.77.gtf";
    let tree = build_tree(path).unwrap();
    let query = Contig::new("FlyBase".to_owned(), 1618414, 4, ReqStrand::Forward);
    tree.find(&query)
        .map(|e| println!("{}", serde_json::to_string_pretty(&e.data()).unwrap()))
        .for_each(drop);
}
