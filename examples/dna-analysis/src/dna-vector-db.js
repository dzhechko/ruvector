/**
 * DNA Vector Database - Ruvector integration for DNA sequence storage and search
 * Provides similarity search and clustering capabilities for genomic data
 */

import { DNAEncoder } from './dna-encoder.js';

/**
 * In-memory vector database implementation for DNA analysis
 * (Simulates ruvector API for pure JavaScript execution)
 */
class VectorDB {
  constructor(options = {}) {
    this.dimensions = options.dimensions || 256;
    this.distanceMetric = options.distanceMetric || 'Cosine';
    this.vectors = new Map();
    this.idCounter = 0;
  }

  async insert(entry) {
    const id = entry.id || `vec_${++this.idCounter}`;
    this.vectors.set(id, {
      id,
      vector: entry.vector,
      metadata: entry.metadata || {}
    });
    return id;
  }

  async insertBatch(entries) {
    const ids = [];
    for (const entry of entries) {
      ids.push(await this.insert(entry));
    }
    return ids;
  }

  async search(query) {
    const results = [];

    for (const [id, entry] of this.vectors) {
      const distance = this.calculateDistance(query.vector, entry.vector);
      results.push({
        id,
        score: distance,
        metadata: entry.metadata,
        vector: query.includeVectors ? entry.vector : undefined
      });
    }

    return results
      .sort((a, b) => a.score - b.score)
      .slice(0, query.k || 10);
  }

  calculateDistance(a, b) {
    if (this.distanceMetric === 'Cosine') {
      return this.cosineDistance(a, b);
    } else if (this.distanceMetric === 'Euclidean') {
      return this.euclideanDistance(a, b);
    }
    return this.cosineDistance(a, b);
  }

  cosineDistance(a, b) {
    let dotProduct = 0;
    let normA = 0;
    let normB = 0;

    for (let i = 0; i < a.length; i++) {
      dotProduct += a[i] * b[i];
      normA += a[i] * a[i];
      normB += b[i] * b[i];
    }

    const similarity = dotProduct / (Math.sqrt(normA) * Math.sqrt(normB) + 1e-10);
    return 1 - similarity; // Convert to distance
  }

  euclideanDistance(a, b) {
    let sum = 0;
    for (let i = 0; i < a.length; i++) {
      const diff = a[i] - b[i];
      sum += diff * diff;
    }
    return Math.sqrt(sum);
  }

  async len() {
    return this.vectors.size;
  }

  async get(id) {
    return this.vectors.get(id) || null;
  }

  async delete(id) {
    return this.vectors.delete(id);
  }

  getAllEntries() {
    return Array.from(this.vectors.values());
  }
}

/**
 * DNA-specific vector database wrapper
 */
export class DNAVectorDB {
  constructor(options = {}) {
    this.kmerSize = options.kmerSize || 4;
    this.encoder = new DNAEncoder(this.kmerSize);

    this.db = new VectorDB({
      dimensions: this.encoder.dimensions,
      distanceMetric: options.distanceMetric || 'Cosine'
    });

    this.sequences = new Map(); // Store original sequences
  }

  /**
   * Store a DNA sequence in the vector database
   */
  async storeSequence(entry) {
    const vector = this.encoder.sequenceToVector(entry.sequence);
    const gcContent = this.encoder.calculateGCContent(entry.sequence);
    const complexity = this.encoder.calculateComplexity(entry.sequence);
    const composition = this.encoder.getNucleotideComposition(entry.sequence);

    const metadata = {
      id: entry.id,
      name: entry.name || entry.id,
      organism: entry.organism || 'Unknown',
      type: entry.type || 'Unknown',
      length: entry.sequence.length,
      gcContent: gcContent.toFixed(2),
      complexity: complexity.toFixed(2),
      composition,
      ...(entry.metadata || {})
    };

    const id = await this.db.insert({
      id: entry.id,
      vector,
      metadata
    });

    this.sequences.set(id, entry.sequence);

    return {
      id,
      metadata,
      vectorDimensions: this.encoder.dimensions
    };
  }

  /**
   * Store multiple sequences in batch
   */
  async storeBatch(sequences) {
    const results = [];
    for (const seq of sequences) {
      results.push(await this.storeSequence(seq));
    }
    return results;
  }

  /**
   * Find similar sequences to a query sequence
   */
  async findSimilar(querySequence, k = 10) {
    const vector = this.encoder.sequenceToVector(querySequence);

    const results = await this.db.search({
      vector,
      k,
      includeVectors: false
    });

    return results.map(r => ({
      ...r,
      similarity: ((1 - r.score) * 100).toFixed(2) + '%'
    }));
  }

  /**
   * Find similar sequences by ID
   */
  async findSimilarById(id, k = 10) {
    const sequence = this.sequences.get(id);
    if (!sequence) {
      throw new Error(`Sequence ${id} not found`);
    }
    return this.findSimilar(sequence, k + 1)
      .then(results => results.filter(r => r.id !== id).slice(0, k));
  }

  /**
   * Compute pairwise similarity matrix
   */
  async computeSimilarityMatrix() {
    const entries = this.db.getAllEntries();
    const n = entries.length;
    const matrix = [];
    const ids = entries.map(e => e.id);

    for (let i = 0; i < n; i++) {
      const row = [];
      for (let j = 0; j < n; j++) {
        if (i === j) {
          row.push(1.0);
        } else {
          const distance = this.db.cosineDistance(
            entries[i].vector,
            entries[j].vector
          );
          row.push(1 - distance);
        }
      }
      matrix.push(row);
    }

    return { ids, matrix };
  }

  /**
   * Cluster sequences using k-means-like algorithm
   */
  async clusterSequences(numClusters = 5, maxIterations = 100) {
    const entries = this.db.getAllEntries();

    if (entries.length < numClusters) {
      numClusters = entries.length;
    }

    // Initialize centroids randomly
    const shuffled = entries.slice().sort(() => Math.random() - 0.5);
    let centroids = shuffled.slice(0, numClusters).map(e =>
      Array.from(e.vector)
    );

    let assignments = new Array(entries.length).fill(0);
    let changed = true;
    let iteration = 0;

    while (changed && iteration < maxIterations) {
      changed = false;
      iteration++;

      // Assign points to nearest centroid
      for (let i = 0; i < entries.length; i++) {
        let minDist = Infinity;
        let bestCluster = 0;

        for (let c = 0; c < numClusters; c++) {
          const dist = this.db.cosineDistance(
            entries[i].vector,
            new Float32Array(centroids[c])
          );
          if (dist < minDist) {
            minDist = dist;
            bestCluster = c;
          }
        }

        if (assignments[i] !== bestCluster) {
          assignments[i] = bestCluster;
          changed = true;
        }
      }

      // Update centroids
      const newCentroids = centroids.map(() =>
        new Array(this.encoder.dimensions).fill(0)
      );
      const counts = new Array(numClusters).fill(0);

      for (let i = 0; i < entries.length; i++) {
        const c = assignments[i];
        counts[c]++;
        for (let d = 0; d < this.encoder.dimensions; d++) {
          newCentroids[c][d] += entries[i].vector[d];
        }
      }

      for (let c = 0; c < numClusters; c++) {
        if (counts[c] > 0) {
          for (let d = 0; d < this.encoder.dimensions; d++) {
            newCentroids[c][d] /= counts[c];
          }
        }
      }

      centroids = newCentroids;
    }

    // Build cluster results
    const clusters = [];
    for (let c = 0; c < numClusters; c++) {
      const members = entries
        .filter((_, i) => assignments[i] === c)
        .map(e => ({
          id: e.id,
          name: e.metadata.name,
          organism: e.metadata.organism,
          gcContent: e.metadata.gcContent
        }));

      if (members.length > 0) {
        clusters.push({
          id: c,
          size: members.length,
          members
        });
      }
    }

    return {
      numClusters: clusters.length,
      iterations: iteration,
      clusters: clusters.sort((a, b) => b.size - a.size)
    };
  }

  /**
   * Analyze a sequence against the database
   */
  async analyzeSequence(sequence, name = 'Query') {
    const gcContent = this.encoder.calculateGCContent(sequence);
    const complexity = this.encoder.calculateComplexity(sequence);
    const composition = this.encoder.getNucleotideComposition(sequence);
    const kmerDist = this.encoder.getKmerDistribution(sequence).slice(0, 20);
    const motifs = this.encoder.findMotifs(sequence);
    const orfs = this.encoder.findORFs(sequence);
    const similar = await this.findSimilar(sequence, 5);

    return {
      name,
      length: sequence.length,
      gcContent: `${gcContent.toFixed(2)}%`,
      complexity: `${complexity.toFixed(2)}%`,
      composition,
      topKmers: kmerDist.map(([kmer, count]) => ({ kmer, count })),
      repeatedMotifs: motifs.map(([motif, count]) => ({ motif, count })),
      openReadingFrames: orfs.slice(0, 5),
      similarSequences: similar
    };
  }

  /**
   * Get database statistics
   */
  async getStats() {
    const entries = this.db.getAllEntries();

    const stats = {
      totalSequences: entries.length,
      totalBasePairs: 0,
      avgLength: 0,
      avgGCContent: 0,
      byOrganism: {},
      byType: {}
    };

    for (const entry of entries) {
      stats.totalBasePairs += entry.metadata.length;
      stats.avgGCContent += parseFloat(entry.metadata.gcContent);

      const org = entry.metadata.organism;
      stats.byOrganism[org] = (stats.byOrganism[org] || 0) + 1;

      const type = entry.metadata.type;
      stats.byType[type] = (stats.byType[type] || 0) + 1;
    }

    if (entries.length > 0) {
      stats.avgLength = Math.round(stats.totalBasePairs / entries.length);
      stats.avgGCContent = (stats.avgGCContent / entries.length).toFixed(2) + '%';
    }

    return stats;
  }
}

export default DNAVectorDB;
