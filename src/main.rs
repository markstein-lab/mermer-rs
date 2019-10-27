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

// The original C implementation of readfasta output several files:
//
// | name           | purpose                                                                       |
// |----------------+-------------------------------------------------------------------------------|
// | genome.txt     | The preprocessed genome (binary).                                             |
// | master.txt     | List of the chromosomes read.                                                 |
// | exceptions.txt | Collection of "exception words". Implementation detail for wildcard searches. |
// | xcontigs.txt   | Chromosome names, indices, locations, and sizes.                              |
// | datasize.txt   | Statistics about the readfasta run.                                           |
//
// The goal of this rewrite is to create a library for the algorithm rather than
// a command-line tool. As such, none of these files are written to disk. The
// data that would be contained in genome.txt, exceptions.txt, and xcontigs.txt
// are all returned as their respective in-memory structures.

use std::error::Error;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;

/// A distinct chromosome in the genome.
#[derive(Debug)]
struct Chromosome {
    name: String,
    position: u64,
    size: u64,
}

/// Buffer used in preprocessing the input .fasta file.
#[derive(Debug)]
struct ParseBuffer {
    nibble: u8,
    nibble_size: u32,
    exception: u32,
    exception_size: u32,

    nibble_backlog: Option<u8>,
    exception_backlog: Option<u32>,
}

impl ParseBuffer {
    /// Return a new, empty ParseBuffer.
    fn new() -> Self {
        ParseBuffer {
            nibble : 0,
            nibble_size: 0,
            exception: 0,
            exception_size: 0,
            nibble_backlog: None,
            exception_backlog: None,
        }
    }

    /// Add a nucleotide to the parse buffer and update the state accordingly.
    ///
    /// # Arguments
    ///
    /// * `nucleotide` - The ASCII code of the nucleotide abbreviation.
    fn add(&mut self, nucleotide: u8) {
        // The mask is used to normalize the case of `nucleotide`.
        match nucleotide & 0xDF {
            0x41 => {
                self.nibble <<= 2;
                self.exception <<= 1;
            }
            0x43 => {
                self.nibble = 1 + (self.nibble << 2);
                self.exception <<= 1;
            }
            0x47 => {
                self.nibble = 2 + (self.nibble << 2);
                self.exception <<= 1;
            }
            0x54 => {
                self.nibble = 3 + (self.nibble << 2);
                self.exception <<= 1;
            }
            _ => {
                // Encode a non-nucleotide character as an 'A' and set a bit in
                // the exception vector.
                self.nibble <<= 2;
                self.exception = 1 + (self.exception << 1);
            }
        }

        self.nibble_size += 1;
        if self.nibble_size == 4 {
            self.nibble_backlog = Some(self.nibble);
            self.nibble = 0;
            self.nibble_size = 0;
        }

        self.exception_size += 1;
        if self.exception_size == 32 {
            if self.exception != 0 {
                self.exception_backlog = Some(self.exception);
            }
            self.exception = 0;
            self.exception_size = 0;
        }
    }

    /// Return the genome nibble if it has been filled.
    fn parsed_nibble(&mut self) -> Option<u8> {
        if let Some(nibble) = self.nibble_backlog {
            self.nibble_backlog = None;
            Some(nibble)
        } else {
            None
        }
    }

    /// Return the exception vector if it has been filled.
    fn parsed_exception(&mut self) -> Option<u32> {
        if let Some(exception) = self.exception_backlog {
            self.exception_backlog = None;
            Some(exception)
        } else {
            None
        }
    }
}

/// Return the preprocessed genome, a list of exceptions, and a list of
/// chromosomes for the genome specified in a .fasta file.
///
/// # Arguments
///
/// * `f` - The .fasta file to read.
///
/// # Example
///
/// ```
/// let f = File::open("dm6.fa").unwrap();
/// let (_, _, chromosomes) = read_fasta(&f).unwrap();
/// for chromosome in chromosomes {
///     println!("{} has {} nucleotides", chromosome.name, chromosome.size);
/// }
/// ```
fn read_fasta(f: &File) -> io::Result<(Vec<u8>, Vec<(u64, u32)>, Vec<Chromosome>)> {
    let mut genome: Vec<u8> = Vec::new();
    let mut exceptions: Vec<(u64, u32)> = Vec::new(); // location, exception
    let mut chromosomes: Vec<Chromosome> = Vec::new();

    let reader = BufReader::new(f);

    let mut position = 0;
    let mut last_position = 0;
    let mut name: Option<String> = None;
    let mut buffer = ParseBuffer::new();

    for line in reader.lines() {
        let line = line?;
        if line.contains(">") {
            if let Some(name) = name {
                chromosomes.push(Chromosome {
                    name,
                    position: last_position,
                    size: position - last_position,
                });
            }

            let index = line.find(">").unwrap() + 1;
            name = Some(String::from(&line[index..]));
            last_position = position;

            // Begin the chromosome with an 'N'.
            buffer.add('N' as u8);
            position += 1;
        } else {
            for nucleotide in line.chars() {
                buffer.add(nucleotide as u8);
                position += 1;
                if let Some(nibble) = buffer.parsed_nibble() {
                    genome.push(nibble);
                }
                if let Some(exception) = buffer.parsed_exception() {
                    exceptions.push(((position - 1) / 32, exception));
                }
            }
        }
    }

    if let Some(name) = name {
        chromosomes.push(Chromosome {
            name,
            position: last_position,
            size: position - last_position,
        });
    }

    // End current chromosome with a N.
    buffer.add('N' as u8);

    loop {
        if let Some(nibble) = buffer.parsed_nibble() {
            genome.push(nibble);
        }
        if let Some(exception) = buffer.parsed_exception() {
            exceptions.push(((position - 1) / 32, exception));
            break;
        }
        buffer.add('N' as u8);
    }

    Ok((genome, exceptions, chromosomes))
}

fn main() {
    let f = File::open("/home/jakob/University/BIOL 296/dm6/dm6.fa").unwrap();
    let (genome, exceptions, chromosomes) = read_fasta(&f).unwrap();
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

    let mut file = match File::create("/tmp/exceptions.txt") {
        Err(why) => panic!("couldn't create {}: {}", "/tmp/test.txt", why.description()),
        Ok(file) => file,
    };

    for (position, word) in exceptions {
        writeln!(file, "{} {:08x}", position, word);
    }
}
