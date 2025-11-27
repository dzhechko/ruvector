/**
 * DNA Data Fetcher - Retrieves public DNA sequences from NCBI GenBank
 */

import fetch from 'node-fetch';
import { writeFileSync, existsSync, readFileSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// NCBI E-utilities base URL
const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

/**
 * Collection of interesting public DNA sequences for analysis
 * Includes various organisms and sequence types
 */
export const PUBLIC_DNA_SEQUENCES = [
  // Human mitochondrial DNA
  { id: 'NC_012920.1', name: 'Homo sapiens mitochondrion', organism: 'Human', type: 'Mitochondrial genome' },

  // Model organism mitochondria
  { id: 'NC_005089.1', name: 'Mus musculus mitochondrion', organism: 'Mouse', type: 'Mitochondrial genome' },
  { id: 'NC_001566.1', name: 'Drosophila melanogaster mitochondrion', organism: 'Fruit fly', type: 'Mitochondrial genome' },
  { id: 'NC_001328.1', name: 'Caenorhabditis elegans mitochondrion', organism: 'Roundworm', type: 'Mitochondrial genome' },

  // Bacterial genomes (partial/small)
  { id: 'NC_000913.3', name: 'Escherichia coli K-12', organism: 'E. coli', type: 'Bacterial chromosome', partial: true },

  // Viral genomes
  { id: 'NC_001802.1', name: 'Human immunodeficiency virus 1', organism: 'HIV-1', type: 'Viral genome' },
  { id: 'NC_045512.2', name: 'SARS-CoV-2 (Wuhan-Hu-1)', organism: 'SARS-CoV-2', type: 'Viral genome' },
  { id: 'NC_006273.2', name: 'Human herpesvirus 5', organism: 'CMV', type: 'Viral genome' },

  // Chloroplast genomes
  { id: 'NC_000932.1', name: 'Arabidopsis thaliana chloroplast', organism: 'Thale cress', type: 'Chloroplast genome' },

  // Famous genes
  { id: 'NM_000546.6', name: 'TP53 tumor protein p53', organism: 'Human', type: 'mRNA' },
  { id: 'NM_005228.5', name: 'EGFR epidermal growth factor receptor', organism: 'Human', type: 'mRNA' },
  { id: 'NM_000518.5', name: 'HBB hemoglobin subunit beta', organism: 'Human', type: 'mRNA' },

  // Additional interesting sequences
  { id: 'NC_001224.1', name: 'Saccharomyces cerevisiae mitochondrion', organism: 'Yeast', type: 'Mitochondrial genome' },
  { id: 'NC_003279.8', name: 'Caenorhabditis elegans chromosome I', organism: 'C. elegans', type: 'Chromosome', partial: true },
];

/**
 * Fetch sequence from NCBI using E-utilities
 */
async function fetchSequenceFromNCBI(accession, retries = 3) {
  const url = `${EUTILS_BASE}/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text`;

  for (let attempt = 0; attempt < retries; attempt++) {
    try {
      console.log(`  Fetching ${accession} (attempt ${attempt + 1})...`);

      const response = await fetch(url, {
        headers: {
          'User-Agent': 'RuvectorDNAAnalysis/1.0 (contact@example.com)'
        }
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const text = await response.text();

      // Parse FASTA format
      const lines = text.trim().split('\n');
      const header = lines[0];
      const sequence = lines.slice(1).join('').replace(/\s/g, '');

      if (sequence.length < 100) {
        throw new Error('Sequence too short or invalid');
      }

      // For very large sequences, take a representative portion
      const maxLength = 100000; // 100kb max
      const finalSequence = sequence.length > maxLength
        ? sequence.substring(0, maxLength)
        : sequence;

      return {
        header: header.substring(1), // Remove '>'
        sequence: finalSequence,
        fullLength: sequence.length,
        truncated: sequence.length > maxLength
      };

    } catch (error) {
      console.log(`    Error: ${error.message}`);
      if (attempt < retries - 1) {
        await new Promise(r => setTimeout(r, 2000 * (attempt + 1)));
      }
    }
  }

  return null;
}

/**
 * Generate synthetic DNA sequences for testing
 */
function generateSyntheticSequences() {
  const nucleotides = ['A', 'C', 'G', 'T'];
  const sequences = [];

  // Random sequences with varying GC content
  const gcLevels = [0.3, 0.4, 0.5, 0.6, 0.7];

  gcLevels.forEach((gcTarget, idx) => {
    let seq = '';
    for (let i = 0; i < 5000; i++) {
      const isGC = Math.random() < gcTarget;
      if (isGC) {
        seq += Math.random() < 0.5 ? 'G' : 'C';
      } else {
        seq += Math.random() < 0.5 ? 'A' : 'T';
      }
    }
    sequences.push({
      id: `SYNTH_GC${Math.round(gcTarget * 100)}`,
      name: `Synthetic sequence GC=${Math.round(gcTarget * 100)}%`,
      organism: 'Synthetic',
      type: 'Synthetic',
      sequence: seq
    });
  });

  // Repetitive sequence
  let repetitive = '';
  const repeat = 'ATCGATCG';
  for (let i = 0; i < 500; i++) {
    repetitive += repeat;
  }
  sequences.push({
    id: 'SYNTH_REPEAT',
    name: 'Synthetic repetitive sequence',
    organism: 'Synthetic',
    type: 'Synthetic',
    sequence: repetitive
  });

  // Coding-like sequence (with start/stop codons)
  let coding = 'ATG'; // Start codon
  for (let i = 0; i < 1000; i++) {
    // Avoid stop codons
    let codon;
    do {
      codon = nucleotides[Math.floor(Math.random() * 4)] +
              nucleotides[Math.floor(Math.random() * 4)] +
              nucleotides[Math.floor(Math.random() * 4)];
    } while (['TAA', 'TAG', 'TGA'].includes(codon));
    coding += codon;
  }
  coding += 'TAA'; // Stop codon
  sequences.push({
    id: 'SYNTH_CODING',
    name: 'Synthetic coding sequence',
    organism: 'Synthetic',
    type: 'Synthetic ORF',
    sequence: coding
  });

  return sequences;
}

/**
 * Fetch all public DNA data
 */
export async function fetchAllDNAData(useCache = true) {
  const cacheFile = join(__dirname, '../data/dna_sequences.json');

  // Check cache
  if (useCache && existsSync(cacheFile)) {
    console.log('Loading cached DNA sequences...');
    const cached = JSON.parse(readFileSync(cacheFile, 'utf-8'));
    console.log(`Loaded ${cached.sequences.length} sequences from cache`);
    return cached;
  }

  console.log('Fetching public DNA sequences from NCBI...\n');

  const sequences = [];

  // Fetch from NCBI
  for (const entry of PUBLIC_DNA_SEQUENCES) {
    const result = await fetchSequenceFromNCBI(entry.id);

    if (result) {
      sequences.push({
        id: entry.id,
        name: entry.name,
        organism: entry.organism,
        type: entry.type,
        header: result.header,
        sequence: result.sequence,
        length: result.sequence.length,
        fullLength: result.fullLength,
        truncated: result.truncated
      });
      console.log(`  ✓ ${entry.name} - ${result.sequence.length} bp`);
    } else {
      console.log(`  ✗ Failed to fetch ${entry.name}`);
    }

    // Rate limiting
    await new Promise(r => setTimeout(r, 500));
  }

  // Add synthetic sequences
  console.log('\nGenerating synthetic sequences...');
  const synthetic = generateSyntheticSequences();
  synthetic.forEach(s => {
    sequences.push({
      ...s,
      length: s.sequence.length,
      fullLength: s.sequence.length,
      truncated: false
    });
    console.log(`  ✓ ${s.name} - ${s.sequence.length} bp`);
  });

  const data = {
    fetchedAt: new Date().toISOString(),
    source: 'NCBI GenBank + Synthetic',
    sequences: sequences
  };

  // Save to cache
  writeFileSync(cacheFile, JSON.stringify(data, null, 2));
  console.log(`\nSaved ${sequences.length} sequences to ${cacheFile}`);

  return data;
}

// Run if executed directly
if (import.meta.url === `file://${process.argv[1]}`) {
  fetchAllDNAData(false).then(data => {
    console.log(`\n✓ Fetched ${data.sequences.length} DNA sequences`);
  }).catch(err => {
    console.error('Error:', err);
    process.exit(1);
  });
}

export default fetchAllDNAData;
