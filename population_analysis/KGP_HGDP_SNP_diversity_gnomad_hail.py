"""Export population-specific AC and AN values for a configurable genomic region."""

import hail as hl # type: ignore
import pandas as pd # type: ignore


HAIL_TABLE = (
    "/scratch/hain/pop_freq/"
    "gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"
)

def main():
    
    
    regions = pd.read_csv("/g/korbel/hain/centromere/revision/HOROSCOPE_manuscript/population_analysis/gr.bed", header=None, names=["chrom", "start", "end"], sep="\t")
    
    for _, row in regions.iterrows():
        REGION = f"{row['chrom']}:{row['start']}-{row['end']}"
        OUTPUT_PATH = f"/scratch/hain/pop_freq/{row['chrom']}_population_ac_an.tsv"
        
        print(row["chrom"])

        try:
            hl.init()

            ht = hl.read_table(HAIL_TABLE)
            region_table = hl.filter_intervals(
                ht,
                [hl.parse_locus_interval(REGION, reference_genome="GRCh38")],
            ).key_by()

            # Keep one adjusted, non-sex-stratified frequency entry per HGDP/1KGP population.
            globals = hl.eval(ht.index_globals())
            population_frequency_indices = {
                metadata["pop"]: index
                for index, metadata in enumerate(globals.hgdp_tgp_freq_meta)
                if metadata.get("group") == "adj"
                and "pop" in metadata
                and "sex" not in metadata
            }

            region_table = region_table.select(
                CHROM=region_table.locus.contig,
                POS=region_table.locus.position,
                REF=region_table.alleles[0],
                ALT=region_table.alleles[1],
                population_counts=hl.array([
                    hl.struct(
                        POPULATION=population,
                        AC=region_table.hgdp_tgp_freq[index].AC,
                        AN=region_table.hgdp_tgp_freq[index].AN,
                    )
                    for population, index in sorted(population_frequency_indices.items())
                ]),
            ).explode("population_counts")
            region_table = region_table.select(
                "CHROM",
                "POS",
                "REF",
                "ALT",
                POPULATION=region_table.population_counts.POPULATION,
                AC=region_table.population_counts.AC,
                AN=region_table.population_counts.AN,
            )

            region_table.export(OUTPUT_PATH)
            hl.stop()
        except Exception as e:
            print(f"Error processing region {REGION}: {e}")
            hl.stop()

if __name__ == "__main__":
    main()