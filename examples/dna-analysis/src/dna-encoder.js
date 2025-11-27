/**
 * DNA Sequence Encoder - Converts DNA sequences to vector embeddings
 * Uses k-mer frequency encoding for vector representation
 */

export class DNAEncoder {
  constructor(kmerSize = 4) {
    this.kmerSize = kmerSize;
    this.kmerMap = this.generateKmerMap();
    this.dimensions = Math.pow(4, kmerSize); // 4^k possible k-mers
  }

  /**
   * Generate mapping of all possible k-mers to indices
   */
  generateKmerMap() {
    const nucleotides = ['A', 'C', 'G', 'T'];
    const kmers = this.generateAllKmers(nucleotides, this.kmerSize);
    const map = new Map();
    kmers.forEach((kmer, index) => map.set(kmer, index));
    return map;
  }

  /**
   * Generate all possible k-mers recursively
   */
  generateAllKmers(nucleotides, k, prefix = '') {
    if (k === 0) return [prefix];
    const result = [];
    for (const n of nucleotides) {
      result.push(...this.generateAllKmers(nucleotides, k - 1, prefix + n));
    }
    return result;
  }

  /**
   * Extract k-mers from a DNA sequence
   */
  extractKmers(sequence) {
    const kmers = [];
    const cleanSeq = sequence.toUpperCase().replace(/[^ACGT]/g, '');

    for (let i = 0; i <= cleanSeq.length - this.kmerSize; i++) {
      const kmer = cleanSeq.substring(i, i + this.kmerSize);
      if (this.kmerMap.has(kmer)) {
        kmers.push(kmer);
      }
    }
    return kmers;
  }

  /**
   * Convert DNA sequence to k-mer frequency vector
   */
  sequenceToVector(sequence) {
    const vector = new Float32Array(this.dimensions);
    const kmers = this.extractKmers(sequence);

    if (kmers.length === 0) {
      return vector;
    }

    // Count k-mer frequencies
    for (const kmer of kmers) {
      const index = this.kmerMap.get(kmer);
      if (index !== undefined) {
        vector[index]++;
      }
    }

    // Normalize to frequencies
    const total = kmers.length;
    for (let i = 0; i < vector.length; i++) {
      vector[i] /= total;
    }

    return vector;
  }

  /**
   * Calculate GC content of a sequence
   */
  calculateGCContent(sequence) {
    const cleanSeq = sequence.toUpperCase().replace(/[^ACGT]/g, '');
    if (cleanSeq.length === 0) return 0;

    const gc = (cleanSeq.match(/[GC]/g) || []).length;
    return (gc / cleanSeq.length) * 100;
  }

  /**
   * Calculate sequence complexity using k-mer diversity
   */
  calculateComplexity(sequence) {
    const kmers = this.extractKmers(sequence);
    const uniqueKmers = new Set(kmers);
    const maxPossible = Math.min(kmers.length, this.dimensions);

    if (maxPossible === 0) return 0;
    return (uniqueKmers.size / maxPossible) * 100;
  }

  /**
   * Get k-mer frequency distribution
   */
  getKmerDistribution(sequence) {
    const kmers = this.extractKmers(sequence);
    const distribution = new Map();

    for (const kmer of kmers) {
      distribution.set(kmer, (distribution.get(kmer) || 0) + 1);
    }

    // Sort by frequency
    return Array.from(distribution.entries())
      .sort((a, b) => b[1] - a[1]);
  }

  /**
   * Find repeated motifs in sequence
   */
  findMotifs(sequence, minLength = 6, minOccurrences = 3) {
    const cleanSeq = sequence.toUpperCase().replace(/[^ACGT]/g, '');
    const motifs = new Map();

    for (let len = minLength; len <= Math.min(20, cleanSeq.length / 2); len++) {
      for (let i = 0; i <= cleanSeq.length - len; i++) {
        const motif = cleanSeq.substring(i, i + len);
        motifs.set(motif, (motifs.get(motif) || 0) + 1);
      }
    }

    return Array.from(motifs.entries())
      .filter(([_, count]) => count >= minOccurrences)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 20);
  }

  /**
   * Calculate nucleotide composition
   */
  getNucleotideComposition(sequence) {
    const cleanSeq = sequence.toUpperCase().replace(/[^ACGT]/g, '');
    const total = cleanSeq.length;

    if (total === 0) return { A: 0, C: 0, G: 0, T: 0 };

    return {
      A: ((cleanSeq.match(/A/g) || []).length / total * 100).toFixed(2),
      C: ((cleanSeq.match(/C/g) || []).length / total * 100).toFixed(2),
      G: ((cleanSeq.match(/G/g) || []).length / total * 100).toFixed(2),
      T: ((cleanSeq.match(/T/g) || []).length / total * 100).toFixed(2)
    };
  }

  /**
   * Detect potential coding regions (simple ORF finder)
   */
  findORFs(sequence, minLength = 100) {
    const cleanSeq = sequence.toUpperCase().replace(/[^ACGT]/g, '');
    const startCodon = 'ATG';
    const stopCodons = ['TAA', 'TAG', 'TGA'];
    const orfs = [];

    for (let frame = 0; frame < 3; frame++) {
      let inORF = false;
      let orfStart = 0;

      for (let i = frame; i < cleanSeq.length - 2; i += 3) {
        const codon = cleanSeq.substring(i, i + 3);

        if (!inORF && codon === startCodon) {
          inORF = true;
          orfStart = i;
        } else if (inORF && stopCodons.includes(codon)) {
          const orfLength = i + 3 - orfStart;
          if (orfLength >= minLength) {
            orfs.push({
              start: orfStart,
              end: i + 3,
              length: orfLength,
              frame: frame + 1,
              sequence: cleanSeq.substring(orfStart, i + 3)
            });
          }
          inORF = false;
        }
      }
    }

    return orfs.sort((a, b) => b.length - a.length);
  }

  /**
   * Encode batch of sequences
   */
  encodeBatch(sequences) {
    return sequences.map(seq => ({
      ...seq,
      vector: this.sequenceToVector(seq.sequence),
      gcContent: this.calculateGCContent(seq.sequence),
      complexity: this.calculateComplexity(seq.sequence),
      composition: this.getNucleotideComposition(seq.sequence)
    }));
  }
}

export default DNAEncoder;
