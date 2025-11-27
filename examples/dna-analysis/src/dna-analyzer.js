/**
 * DNA Analyzer - Comprehensive DNA sequence analysis module
 * Provides phylogenetic analysis, mutation detection, and visualization
 */

import { DNAEncoder } from './dna-encoder.js';

export class DNAAnalyzer {
  constructor() {
    this.encoder = new DNAEncoder(4);
  }

  /**
   * Calculate Levenshtein edit distance between two sequences
   */
  calculateEditDistance(seq1, seq2, maxLength = 5000) {
    // Truncate for performance
    const s1 = seq1.substring(0, maxLength);
    const s2 = seq2.substring(0, maxLength);

    const m = s1.length;
    const n = s2.length;

    // Use space-efficient version
    let prev = new Array(n + 1).fill(0).map((_, i) => i);
    let curr = new Array(n + 1).fill(0);

    for (let i = 1; i <= m; i++) {
      curr[0] = i;
      for (let j = 1; j <= n; j++) {
        if (s1[i - 1] === s2[j - 1]) {
          curr[j] = prev[j - 1];
        } else {
          curr[j] = 1 + Math.min(prev[j - 1], prev[j], curr[j - 1]);
        }
      }
      [prev, curr] = [curr, prev];
    }

    return prev[n];
  }

  /**
   * Find mutations between two aligned sequences
   */
  findMutations(seq1, seq2, windowSize = 100) {
    const mutations = [];
    const minLen = Math.min(seq1.length, seq2.length);

    for (let i = 0; i < minLen; i++) {
      if (seq1[i] !== seq2[i]) {
        mutations.push({
          position: i,
          original: seq1[i],
          mutated: seq2[i],
          type: this.classifyMutation(seq1[i], seq2[i])
        });
      }
    }

    // Summarize by region
    const regions = [];
    for (let i = 0; i < minLen; i += windowSize) {
      const regionMutations = mutations.filter(
        m => m.position >= i && m.position < i + windowSize
      );
      regions.push({
        start: i,
        end: Math.min(i + windowSize, minLen),
        mutationCount: regionMutations.length,
        mutationRate: (regionMutations.length / windowSize * 100).toFixed(2) + '%'
      });
    }

    return {
      totalMutations: mutations.length,
      mutationRate: ((mutations.length / minLen) * 100).toFixed(4) + '%',
      mutations: mutations.slice(0, 100), // Limit for display
      hotspots: regions
        .filter(r => r.mutationCount > 0)
        .sort((a, b) => b.mutationCount - a.mutationCount)
        .slice(0, 10)
    };
  }

  /**
   * Classify mutation type (transition vs transversion)
   */
  classifyMutation(from, to) {
    const purines = ['A', 'G'];
    const pyrimidines = ['C', 'T'];

    const fromPurine = purines.includes(from);
    const toPurine = purines.includes(to);

    if (fromPurine === toPurine) {
      return 'transition'; // Purine↔Purine or Pyrimidine↔Pyrimidine
    }
    return 'transversion'; // Purine↔Pyrimidine
  }

  /**
   * Compute codon usage bias
   */
  analyzeCodonUsage(sequence) {
    const cleanSeq = sequence.toUpperCase().replace(/[^ACGT]/g, '');
    const codonCounts = new Map();

    // Standard genetic code
    const geneticCode = {
      'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
      'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
      'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
      'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
      'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
      'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
      'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
      'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
      'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
      'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
      'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
      'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
      'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
      'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
      'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
      'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    };

    // Count codons
    for (let i = 0; i <= cleanSeq.length - 3; i += 3) {
      const codon = cleanSeq.substring(i, i + 3);
      if (geneticCode[codon]) {
        codonCounts.set(codon, (codonCounts.get(codon) || 0) + 1);
      }
    }

    // Group by amino acid
    const aminoAcidCodons = {};
    for (const [codon, count] of codonCounts) {
      const aa = geneticCode[codon];
      if (!aminoAcidCodons[aa]) {
        aminoAcidCodons[aa] = [];
      }
      aminoAcidCodons[aa].push({ codon, count });
    }

    // Calculate relative synonymous codon usage (RSCU)
    const rscu = {};
    for (const [aa, codons] of Object.entries(aminoAcidCodons)) {
      const total = codons.reduce((sum, c) => sum + c.count, 0);
      const expected = total / codons.length;
      rscu[aa] = codons.map(c => ({
        codon: c.codon,
        count: c.count,
        rscu: expected > 0 ? (c.count / expected).toFixed(3) : 0
      }));
    }

    return {
      totalCodons: Array.from(codonCounts.values()).reduce((a, b) => a + b, 0),
      codonCounts: Object.fromEntries(
        Array.from(codonCounts.entries())
          .sort((a, b) => b[1] - a[1])
      ),
      rscu
    };
  }

  /**
   * Detect CpG islands (methylation-prone regions)
   */
  findCpGIslands(sequence, windowSize = 200, minCpG = 0.6, minGC = 0.5) {
    const cleanSeq = sequence.toUpperCase().replace(/[^ACGT]/g, '');
    const islands = [];

    for (let i = 0; i <= cleanSeq.length - windowSize; i += 50) {
      const window = cleanSeq.substring(i, i + windowSize);

      // Calculate observed CpG
      const cpgCount = (window.match(/CG/g) || []).length;

      // Calculate expected CpG
      const cCount = (window.match(/C/g) || []).length;
      const gCount = (window.match(/G/g) || []).length;
      const expectedCpG = (cCount * gCount) / windowSize;

      // Calculate GC content
      const gcContent = (cCount + gCount) / windowSize;

      // CpG ratio
      const cpgRatio = expectedCpG > 0 ? cpgCount / expectedCpG : 0;

      if (cpgRatio >= minCpG && gcContent >= minGC) {
        // Check if extends previous island
        if (islands.length > 0 && i <= islands[islands.length - 1].end + 50) {
          islands[islands.length - 1].end = i + windowSize;
        } else {
          islands.push({
            start: i,
            end: i + windowSize,
            cpgRatio: cpgRatio.toFixed(3),
            gcContent: (gcContent * 100).toFixed(2) + '%'
          });
        }
      }
    }

    return islands;
  }

  /**
   * Generate phylogenetic distance matrix
   */
  computePhylogeneticDistances(sequences) {
    const n = sequences.length;
    const matrix = [];
    const names = sequences.map(s => s.name || s.id);

    for (let i = 0; i < n; i++) {
      const row = [];
      for (let j = 0; j < n; j++) {
        if (i === j) {
          row.push(0);
        } else if (j < i) {
          row.push(matrix[j][i]); // Symmetric
        } else {
          // Use k-mer distance for speed
          const v1 = this.encoder.sequenceToVector(sequences[i].sequence);
          const v2 = this.encoder.sequenceToVector(sequences[j].sequence);
          const distance = this.cosineDistance(v1, v2);
          row.push(parseFloat(distance.toFixed(4)));
        }
      }
      matrix.push(row);
    }

    return { names, matrix };
  }

  cosineDistance(a, b) {
    let dot = 0, normA = 0, normB = 0;
    for (let i = 0; i < a.length; i++) {
      dot += a[i] * b[i];
      normA += a[i] * a[i];
      normB += b[i] * b[i];
    }
    return 1 - (dot / (Math.sqrt(normA) * Math.sqrt(normB) + 1e-10));
  }

  /**
   * Build simple UPGMA tree
   */
  buildUPGMATree(sequences) {
    const { names, matrix } = this.computePhylogeneticDistances(sequences);
    const n = names.length;

    // Initialize clusters
    const clusters = names.map((name, i) => ({
      name,
      height: 0,
      members: [i],
      isLeaf: true
    }));

    const distances = matrix.map(row => [...row]);
    const active = new Set(Array.from({ length: n }, (_, i) => i));

    let nodeCounter = n;

    while (active.size > 1) {
      // Find minimum distance pair
      let minDist = Infinity;
      let minI = -1, minJ = -1;

      for (const i of active) {
        for (const j of active) {
          if (i < j && distances[i][j] < minDist) {
            minDist = distances[i][j];
            minI = i;
            minJ = j;
          }
        }
      }

      if (minI === -1) break;

      // Merge clusters
      const newCluster = {
        name: `Node${nodeCounter++}`,
        height: minDist / 2,
        children: [clusters[minI], clusters[minJ]],
        members: [...clusters[minI].members, ...clusters[minJ].members],
        isLeaf: false
      };

      // Update distances (UPGMA averaging)
      const newIndex = minI;
      const ni = clusters[minI].members.length;
      const nj = clusters[minJ].members.length;

      for (const k of active) {
        if (k !== minI && k !== minJ) {
          const avgDist = (distances[minI][k] * ni + distances[minJ][k] * nj) / (ni + nj);
          distances[newIndex][k] = avgDist;
          distances[k][newIndex] = avgDist;
        }
      }

      clusters[minI] = newCluster;
      active.delete(minJ);
    }

    // Return root
    const rootIndex = active.values().next().value;
    return clusters[rootIndex];
  }

  /**
   * Convert tree to Newick format
   */
  treeToNewick(node) {
    if (node.isLeaf) {
      return node.name;
    }

    const children = node.children.map(c => {
      const childNewick = this.treeToNewick(c);
      const branchLength = (node.height - c.height).toFixed(4);
      return `${childNewick}:${branchLength}`;
    });

    return `(${children.join(',')})`;
  }

  /**
   * Generate ASCII tree visualization
   */
  treeToASCII(node, prefix = '', isLast = true) {
    let result = prefix;
    result += isLast ? '└── ' : '├── ';
    result += node.isLeaf ? node.name : `[${node.members.length} seqs]`;
    result += '\n';

    if (!node.isLeaf && node.children) {
      const childPrefix = prefix + (isLast ? '    ' : '│   ');
      node.children.forEach((child, i) => {
        result += this.treeToASCII(
          child,
          childPrefix,
          i === node.children.length - 1
        );
      });
    }

    return result;
  }

  /**
   * Perform comprehensive sequence analysis
   */
  analyzeSequence(sequence, name = 'Sequence') {
    const gcContent = this.encoder.calculateGCContent(sequence);
    const complexity = this.encoder.calculateComplexity(sequence);
    const composition = this.encoder.getNucleotideComposition(sequence);
    const kmerDist = this.encoder.getKmerDistribution(sequence);
    const motifs = this.encoder.findMotifs(sequence);
    const orfs = this.encoder.findORFs(sequence);
    const cpgIslands = this.findCpGIslands(sequence);
    const codonUsage = this.analyzeCodonUsage(sequence);

    return {
      name,
      length: sequence.length,
      gcContent: `${gcContent.toFixed(2)}%`,
      complexity: `${complexity.toFixed(2)}%`,
      composition,
      topKmers: kmerDist.slice(0, 15).map(([kmer, count]) => ({ kmer, count })),
      repeatedMotifs: motifs.slice(0, 10),
      openReadingFrames: orfs.slice(0, 5),
      cpgIslands: cpgIslands.slice(0, 10),
      codonUsage: {
        totalCodons: codonUsage.totalCodons,
        topCodons: Object.entries(codonUsage.codonCounts)
          .slice(0, 10)
          .map(([codon, count]) => ({ codon, count }))
      }
    };
  }
}

export default DNAAnalyzer;
