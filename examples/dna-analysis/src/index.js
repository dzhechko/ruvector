#!/usr/bin/env node
/**
 * RuVector DNA Analysis - Main Entry Point
 * Performs comprehensive DNA sequence analysis using vector embeddings
 */

import { fetchAllDNAData, PUBLIC_DNA_SEQUENCES } from './fetch-dna-data.js';
import { DNAVectorDB } from './dna-vector-db.js';
import { DNAAnalyzer } from './dna-analyzer.js';
import { writeFileSync, mkdirSync, existsSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

const RESULTS_DIR = join(__dirname, '../results');

/**
 * Format number with thousands separator
 */
function formatNumber(n) {
  return n.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
}

/**
 * Print section header
 */
function printSection(title) {
  console.log('\n' + '═'.repeat(70));
  console.log(`  ${title}`);
  console.log('═'.repeat(70));
}

/**
 * Print subsection
 */
function printSubsection(title) {
  console.log(`\n--- ${title} ---`);
}

/**
 * Main analysis pipeline
 */
async function runDNAAnalysis() {
  console.log(`
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║        🧬 RuVector DNA Analysis Pipeline 🧬                          ║
║                                                                      ║
║   High-Performance Vector Database for Genomic Analysis              ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
`);

  // Ensure results directory exists
  if (!existsSync(RESULTS_DIR)) {
    mkdirSync(RESULTS_DIR, { recursive: true });
  }

  const analysisResults = {
    timestamp: new Date().toISOString(),
    pipeline: 'RuVector DNA Analysis',
    version: '1.0.0',
    sections: {}
  };

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 1: Fetch DNA Data
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 1: Fetching Public DNA Sequences');

  let dnaData;
  try {
    dnaData = await fetchAllDNAData(true);
    console.log(`\n✓ Loaded ${dnaData.sequences.length} DNA sequences`);
    console.log(`  Source: ${dnaData.source}`);
    console.log(`  Fetched: ${dnaData.fetchedAt}`);
  } catch (error) {
    console.log('Could not fetch from NCBI, generating sample data...');
    // Generate sample data if fetch fails
    dnaData = generateSampleData();
  }

  analysisResults.sections.dataSource = {
    source: dnaData.source,
    fetchedAt: dnaData.fetchedAt,
    sequenceCount: dnaData.sequences.length
  };

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 2: Initialize Vector Database
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 2: Initializing RuVector Database');

  const vectorDB = new DNAVectorDB({
    kmerSize: 4,
    distanceMetric: 'Cosine'
  });

  console.log(`\n✓ Database initialized`);
  console.log(`  K-mer size: 4`);
  console.log(`  Vector dimensions: ${vectorDB.encoder.dimensions} (4^4 = 256)`);
  console.log(`  Distance metric: Cosine similarity`);

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 3: Encode and Store Sequences
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 3: Encoding DNA Sequences to Vectors');

  console.log('\nProcessing sequences...\n');

  const storedSequences = [];
  for (const seq of dnaData.sequences) {
    const result = await vectorDB.storeSequence(seq);
    storedSequences.push(result);

    console.log(`  ✓ ${seq.name.substring(0, 40).padEnd(40)} | ${formatNumber(seq.length).padStart(8)} bp | GC: ${result.metadata.gcContent}%`);
  }

  const dbStats = await vectorDB.getStats();
  console.log(`\n✓ Stored ${dbStats.totalSequences} sequences`);
  console.log(`  Total base pairs: ${formatNumber(dbStats.totalBasePairs)}`);
  console.log(`  Average length: ${formatNumber(dbStats.avgLength)} bp`);
  console.log(`  Average GC content: ${dbStats.avgGCContent}`);

  analysisResults.sections.database = dbStats;

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 4: Similarity Analysis
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 4: K-mer Based Similarity Analysis');

  console.log('\nComputing pairwise similarities...');

  const similarityMatrix = await vectorDB.computeSimilarityMatrix();

  // Find most similar pairs
  const pairs = [];
  for (let i = 0; i < similarityMatrix.ids.length; i++) {
    for (let j = i + 1; j < similarityMatrix.ids.length; j++) {
      pairs.push({
        seq1: similarityMatrix.ids[i],
        seq2: similarityMatrix.ids[j],
        similarity: similarityMatrix.matrix[i][j]
      });
    }
  }

  const topSimilar = pairs.sort((a, b) => b.similarity - a.similarity).slice(0, 10);

  printSubsection('Top 10 Most Similar Sequence Pairs');
  console.log();
  for (const pair of topSimilar) {
    const name1 = (await vectorDB.db.get(pair.seq1))?.metadata?.name?.substring(0, 25) || pair.seq1;
    const name2 = (await vectorDB.db.get(pair.seq2))?.metadata?.name?.substring(0, 25) || pair.seq2;
    console.log(`  ${(pair.similarity * 100).toFixed(2)}% | ${name1.padEnd(25)} <-> ${name2}`);
  }

  analysisResults.sections.similarity = {
    topPairs: await Promise.all(topSimilar.slice(0, 10).map(async p => ({
      sequence1: (await vectorDB.db.get(p.seq1))?.metadata?.name || p.seq1,
      sequence2: (await vectorDB.db.get(p.seq2))?.metadata?.name || p.seq2,
      similarity: `${(p.similarity * 100).toFixed(2)}%`
    })))
  };

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 5: Sequence Clustering
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 5: Sequence Clustering');

  console.log('\nPerforming k-means clustering on vector embeddings...');

  const clusters = await vectorDB.clusterSequences(5);

  console.log(`\n✓ Found ${clusters.numClusters} clusters in ${clusters.iterations} iterations\n`);

  for (const cluster of clusters.clusters) {
    console.log(`  Cluster ${cluster.id + 1} (${cluster.size} sequences):`);
    for (const member of cluster.members.slice(0, 5)) {
      console.log(`    - ${member.name} [${member.organism}]`);
    }
    if (cluster.members.length > 5) {
      console.log(`    ... and ${cluster.members.length - 5} more`);
    }
    console.log();
  }

  analysisResults.sections.clustering = clusters;

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 6: Detailed Sequence Analysis
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 6: Detailed Sequence Analysis');

  const analyzer = new DNAAnalyzer();

  // Analyze key sequences
  const keySequences = dnaData.sequences.filter(s =>
    s.type.includes('Viral') ||
    s.type.includes('Mitochondrial') ||
    s.organism === 'Human'
  ).slice(0, 5);

  const detailedAnalyses = [];

  for (const seq of keySequences) {
    console.log(`\n--- ${seq.name} ---`);

    const analysis = await vectorDB.analyzeSequence(seq.sequence, seq.name);
    detailedAnalyses.push(analysis);

    console.log(`  Length: ${formatNumber(seq.sequence.length)} bp`);
    console.log(`  GC Content: ${analysis.gcContent}`);
    console.log(`  Sequence Complexity: ${analysis.complexity}`);
    console.log(`  Composition: A=${analysis.composition.A}% C=${analysis.composition.C}% G=${analysis.composition.G}% T=${analysis.composition.T}%`);

    if (analysis.openReadingFrames.length > 0) {
      console.log(`  Open Reading Frames: ${analysis.openReadingFrames.length} detected`);
      console.log(`    Longest ORF: ${formatNumber(analysis.openReadingFrames[0].length)} bp (frame ${analysis.openReadingFrames[0].frame})`);
    }

    if (analysis.repeatedMotifs.length > 0) {
      console.log(`  Repeated Motifs: ${analysis.repeatedMotifs.length} found`);
      const topMotif = analysis.repeatedMotifs[0];
      console.log(`    Top motif: "${topMotif[0]}" (${topMotif[1]} occurrences)`);
    }

    console.log(`  Most Similar:`);
    for (const sim of analysis.similarSequences.slice(0, 3)) {
      console.log(`    - ${sim.metadata.name}: ${sim.similarity}`);
    }
  }

  analysisResults.sections.detailedAnalysis = detailedAnalyses;

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 7: Phylogenetic Analysis
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 7: Phylogenetic Analysis');

  // Use mitochondrial sequences for phylogeny
  const mitoSequences = dnaData.sequences.filter(s =>
    s.type.includes('Mitochondrial')
  );

  if (mitoSequences.length >= 3) {
    console.log(`\nBuilding UPGMA tree from ${mitoSequences.length} mitochondrial sequences...\n`);

    const tree = analyzer.buildUPGMATree(mitoSequences);
    const asciiTree = analyzer.treeToASCII(tree);
    const newick = analyzer.treeToNewick(tree) + ';';

    console.log('Phylogenetic Tree (ASCII):');
    console.log(asciiTree);
    console.log('Newick Format:');
    console.log(newick);

    analysisResults.sections.phylogeny = {
      method: 'UPGMA',
      sequenceCount: mitoSequences.length,
      newick,
      asciiTree
    };
  }

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 8: Comparative Genomics
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 8: Comparative Genomics');

  // Compare viral genomes
  const viralSequences = dnaData.sequences.filter(s => s.type.includes('Viral'));

  if (viralSequences.length >= 2) {
    console.log(`\nComparing ${viralSequences.length} viral genomes...\n`);

    const comparisons = [];
    for (let i = 0; i < viralSequences.length; i++) {
      for (let j = i + 1; j < viralSequences.length; j++) {
        const v1 = vectorDB.encoder.sequenceToVector(viralSequences[i].sequence);
        const v2 = vectorDB.encoder.sequenceToVector(viralSequences[j].sequence);

        let dot = 0, n1 = 0, n2 = 0;
        for (let k = 0; k < v1.length; k++) {
          dot += v1[k] * v2[k];
          n1 += v1[k] * v1[k];
          n2 += v2[k] * v2[k];
        }
        const similarity = dot / (Math.sqrt(n1) * Math.sqrt(n2) + 1e-10);

        comparisons.push({
          virus1: viralSequences[i].name,
          virus2: viralSequences[j].name,
          similarity: similarity,
          gc1: vectorDB.encoder.calculateGCContent(viralSequences[i].sequence),
          gc2: vectorDB.encoder.calculateGCContent(viralSequences[j].sequence)
        });
      }
    }

    for (const comp of comparisons) {
      console.log(`  ${comp.virus1.substring(0, 30).padEnd(30)} vs ${comp.virus2.substring(0, 30)}`);
      console.log(`    K-mer similarity: ${(comp.similarity * 100).toFixed(2)}%`);
      console.log(`    GC content: ${comp.gc1.toFixed(1)}% vs ${comp.gc2.toFixed(1)}%`);
      console.log();
    }

    analysisResults.sections.viralComparison = comparisons;
  }

  // ═══════════════════════════════════════════════════════════════════
  // PHASE 9: Summary Statistics
  // ═══════════════════════════════════════════════════════════════════
  printSection('PHASE 9: Summary Statistics');

  console.log('\n--- Database Overview ---');
  console.log(`  Total sequences: ${dbStats.totalSequences}`);
  console.log(`  Total base pairs: ${formatNumber(dbStats.totalBasePairs)}`);
  console.log(`  Average GC content: ${dbStats.avgGCContent}`);

  console.log('\n--- By Organism ---');
  for (const [org, count] of Object.entries(dbStats.byOrganism)) {
    console.log(`  ${org.padEnd(20)}: ${count} sequences`);
  }

  console.log('\n--- By Type ---');
  for (const [type, count] of Object.entries(dbStats.byType)) {
    console.log(`  ${type.padEnd(20)}: ${count} sequences`);
  }

  // ═══════════════════════════════════════════════════════════════════
  // Save Results
  // ═══════════════════════════════════════════════════════════════════
  printSection('Saving Analysis Results');

  const resultsFile = join(RESULTS_DIR, `analysis_${Date.now()}.json`);
  writeFileSync(resultsFile, JSON.stringify(analysisResults, null, 2));
  console.log(`\n✓ Results saved to: ${resultsFile}`);

  // Generate summary report
  const summaryReport = generateSummaryReport(analysisResults, dbStats);
  const summaryFile = join(RESULTS_DIR, `summary_${Date.now()}.txt`);
  writeFileSync(summaryFile, summaryReport);
  console.log(`✓ Summary saved to: ${summaryFile}`);

  console.log(`
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║               🧬 Analysis Complete! 🧬                               ║
║                                                                      ║
║   Results saved in: examples/dna-analysis/results/                   ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
`);

  return analysisResults;
}

/**
 * Generate sample data when NCBI fetch fails
 */
function generateSampleData() {
  const nucleotides = ['A', 'C', 'G', 'T'];

  function randomSeq(length, gcBias = 0.5) {
    let seq = '';
    for (let i = 0; i < length; i++) {
      const isGC = Math.random() < gcBias;
      if (isGC) {
        seq += Math.random() < 0.5 ? 'G' : 'C';
      } else {
        seq += Math.random() < 0.5 ? 'A' : 'T';
      }
    }
    return seq;
  }

  return {
    fetchedAt: new Date().toISOString(),
    source: 'Synthetic (NCBI unavailable)',
    sequences: [
      { id: 'SAMPLE_HUMAN_MITO', name: 'Sample Human Mitochondrion', organism: 'Human', type: 'Mitochondrial genome', sequence: randomSeq(16569, 0.44) },
      { id: 'SAMPLE_MOUSE_MITO', name: 'Sample Mouse Mitochondrion', organism: 'Mouse', type: 'Mitochondrial genome', sequence: randomSeq(16299, 0.37) },
      { id: 'SAMPLE_YEAST_MITO', name: 'Sample Yeast Mitochondrion', organism: 'Yeast', type: 'Mitochondrial genome', sequence: randomSeq(85779, 0.17) },
      { id: 'SAMPLE_SARS2', name: 'Sample SARS-CoV-2', organism: 'SARS-CoV-2', type: 'Viral genome', sequence: randomSeq(29903, 0.38) },
      { id: 'SAMPLE_HIV', name: 'Sample HIV-1', organism: 'HIV-1', type: 'Viral genome', sequence: randomSeq(9719, 0.41) },
      { id: 'SAMPLE_TP53', name: 'Sample TP53 mRNA', organism: 'Human', type: 'mRNA', sequence: randomSeq(2591, 0.47) },
      { id: 'SAMPLE_ECOLI', name: 'Sample E. coli (partial)', organism: 'E. coli', type: 'Bacterial chromosome', sequence: randomSeq(50000, 0.51) },
    ]
  };
}

/**
 * Generate text summary report
 */
function generateSummaryReport(results, stats) {
  let report = `
═══════════════════════════════════════════════════════════════════════════════
                         RUVECTOR DNA ANALYSIS REPORT
═══════════════════════════════════════════════════════════════════════════════

Generated: ${results.timestamp}
Pipeline: ${results.pipeline} v${results.version}

───────────────────────────────────────────────────────────────────────────────
                              DATABASE SUMMARY
───────────────────────────────────────────────────────────────────────────────

Total Sequences: ${stats.totalSequences}
Total Base Pairs: ${formatNumber(stats.totalBasePairs)}
Average Length: ${formatNumber(stats.avgLength)} bp
Average GC Content: ${stats.avgGCContent}

By Organism:
${Object.entries(stats.byOrganism).map(([k, v]) => `  - ${k}: ${v}`).join('\n')}

By Type:
${Object.entries(stats.byType).map(([k, v]) => `  - ${k}: ${v}`).join('\n')}

───────────────────────────────────────────────────────────────────────────────
                            SIMILARITY ANALYSIS
───────────────────────────────────────────────────────────────────────────────

Top Similar Sequence Pairs:
${results.sections.similarity?.topPairs?.map((p, i) =>
  `  ${i + 1}. ${p.sequence1} <-> ${p.sequence2}: ${p.similarity}`
).join('\n') || 'N/A'}

───────────────────────────────────────────────────────────────────────────────
                             CLUSTER ANALYSIS
───────────────────────────────────────────────────────────────────────────────

Number of Clusters: ${results.sections.clustering?.numClusters || 'N/A'}
Iterations: ${results.sections.clustering?.iterations || 'N/A'}

${results.sections.clustering?.clusters?.map(c =>
  `Cluster ${c.id + 1} (${c.size} sequences):\n` +
  c.members.slice(0, 5).map(m => `  - ${m.name} [${m.organism}]`).join('\n')
).join('\n\n') || 'N/A'}

───────────────────────────────────────────────────────────────────────────────
                           PHYLOGENETIC ANALYSIS
───────────────────────────────────────────────────────────────────────────────

Method: ${results.sections.phylogeny?.method || 'N/A'}
Sequences Analyzed: ${results.sections.phylogeny?.sequenceCount || 'N/A'}

${results.sections.phylogeny?.asciiTree || 'Tree visualization not available'}

───────────────────────────────────────────────────────────────────────────────
                                   NOTES
───────────────────────────────────────────────────────────────────────────────

- K-mer size: 4 (256 dimensions)
- Distance metric: Cosine similarity
- Clustering algorithm: K-means
- Phylogeny method: UPGMA (Unweighted Pair Group Method with Arithmetic Mean)

For detailed JSON results, see the accompanying analysis_*.json file.

═══════════════════════════════════════════════════════════════════════════════
                              END OF REPORT
═══════════════════════════════════════════════════════════════════════════════
`;

  return report;
}

// Run analysis
runDNAAnalysis().catch(err => {
  console.error('Analysis failed:', err);
  process.exit(1);
});
