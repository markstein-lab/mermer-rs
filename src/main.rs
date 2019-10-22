// Copyright © 2019 Jakob L. Kreuze <zerodaysfordays@sdf.lonestar.org>
// 
// This file is part of mermer-rs.
// 
// mermer-rs is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
// 
// mermer-rs is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mermer-rs. If not, see <http://www.gnu.org/licenses/>.

// | name           | purpose                                                                       |
// |----------------+-------------------------------------------------------------------------------|
// | genome.txt     | The preprocessed genome (binary).                                             |
// | master.txt     | List of the chromosomes read.                                                 |
// | exceptions.txt | Collection of "exception words". Implementation detail for wildcard searches. |
// | xcontigs.txt   | Chromosome names, indices, locations, and sizes.                              |
// | datasize.txt   | Statistics about the readfasta run.                                           |

use std::error::Error;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;

#[derive(Debug)]
struct Chromosome {
    name: String,
    position: u64,
    size: u64,
}

// TODO: Return exception too.
fn read_fasta(f: &File) -> io::Result<(Vec<u8>, Vec<Chromosome>)> {
    let reader = BufReader::new(f);
    let mut genome: Vec<u8> = Vec::new();
    let mut exceptions: Vec<(u64, u8)> = Vec::new(); // location, exception
    let mut result: Vec<Chromosome> = Vec::new();

    let mut position = 0;
    let mut last_position = 0;
    let mut name: Option<String> = None;

    let mut genome_buffer = 0;
    let mut genome_buffer_size = 0;
    let mut exception_buffer = 0;
    let mut exception_buffer_size = 0;

    for line in reader.lines() {
        let line = line?;
        if line.contains(">") {
            if let Some(name) = name {
                result.push(Chromosome {
                    name,
                    position: last_position,
                    size: position - last_position,
                });
            }

            let index = line.find(">").unwrap() + 1;
            name = Some(String::from(&line[index..]));
            last_position = position;

            // Start by inserting an 'N' into the genome.
            genome_buffer <<= 2;
            genome_buffer_size += 1;
            exception_buffer = 1 + (exception_buffer << 1);
            exception_buffer_size += 1;
            position += 1;
        } else {
            for nucleotide in line.chars() {
                match (nucleotide as u8) & 0xDF {
                    0x41 => {
                        genome_buffer <<= 2;
                        exception_buffer <<= 1;
                    }
                    0x43 => {
                        genome_buffer = 1 + (genome_buffer << 2);
                        exception_buffer <<= 1;
                    }
                    0x47 => {
                        genome_buffer = 2 + (genome_buffer << 2);
                        exception_buffer <<= 1;
                    }
                    0x54 => {
                        genome_buffer = 3 + (genome_buffer << 2);
                        exception_buffer <<= 1;
                    }
                    _ => {
                        genome_buffer <<= 2;
                        exception_buffer = 1 + (exception_buffer << 1);
                    }
                }

                genome_buffer_size += 1;
                if genome_buffer_size == 4 {
                    genome.push(genome_buffer);
                    genome_buffer = 0;
                    genome_buffer_size = 0;
                }

                exception_buffer_size += 1;
                if exception_buffer_size == 32 {
                    if exception_buffer != 0 {
                        exceptions.push((position, exception_buffer));
                    }
                    exception_buffer = 0;
                    exception_buffer_size = 0;
                }

                position += 1;
            }
        }
    }

    if let Some(name) = name {
        result.push(Chromosome {
            name,
            position: last_position,
            size: position - last_position,
        });
    }

    // End current chromosome with a N.
    genome_buffer <<= 2;
    genome_buffer_size += 1;
    if genome_buffer_size == 4 {
        genome.push(genome_buffer);
        genome_buffer = 0;
        genome_buffer_size = 0;
    }

    exception_buffer = 1 + (exception_buffer << 1);
    exception_buffer_size += 1;
    if exception_buffer_size == 32 {
        if exception_buffer != 0 {
            exceptions.push((position, exception_buffer));
        }
        exception_buffer = 0;
        exception_buffer_size = 0;
    }
    
    position += 1;

    while exception_buffer_size != 32 {
        genome_buffer <<= 2;
        exception_buffer = 1 + (exception_buffer << 1);

        genome_buffer_size += 1;
        if genome_buffer_size == 4 {
            genome.push(genome_buffer);
            genome_buffer = 0;
            genome_buffer_size = 0;
        }

        exception_buffer_size += 1;

        position += 1;
    }
    
    if exception_buffer != 0 {
        exceptions.push((position, exception_buffer));
    }

    Ok((genome, result))
}

fn main() {
    let f = File::open("/home/jakob/University/BIOL 296/dm6/dm6.fa").unwrap();
    let (genome, chromosomes) = read_fasta(&f).unwrap();
    for chromosome in chromosomes {
        println!("{:?}", chromosome);
    }

    let mut file = match File::create("/tmp/test.txt") {
        Err(why) => panic!("couldn't create {}: {}", "/tmp/test.txt", why.description()),
        Ok(file) => file,
    };

    match file.write_all(&genome) {
        Err(why) => panic!("couldn't write to {}: {}", "/tmp/test.txt", why.description()),
        Ok(_) => println!("successfully wrote to {}", "/tmp/test.txt"),
    };
}
