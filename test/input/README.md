# test input log

> **Note on generated VCF fixtures**
> `1_labelled_variant.vcf.bgz` and `newphase.vcf.bgz` are **not committed** to the
> repository (binary + index files). They are synthetic and must be generated once:
>
>     python3 test/input/create_test_fixtures.py
>
> Requires `bgzip` and `tabix` (part of htslib). Tests that depend on these files
> will be skipped with an informative message until the fixtures are present.

#### 1_labelled_variant.vcf.bgz  *(generated — not committed)*

- Two variants on chr20 (63406931 C>CGG, 63406991 C>CGG)
- Both annotated against ENSG00000075043
- Variant 1: CategoryBoolean3=1, CategorySample4=male
- Variant 2: CategoryBoolean3=1, CategoryBoolean1=1, CategorySample4=male
- Samples: male (het, 0/1), mother_1 (WT), father_1 (WT)

#### aip_output_example.json

- Single variant, single sample, minimised output
- Most VCF Info fields removed for size
- Sample ID: SAMPLE_1
- variant = 7-93105286-T-A__ENSG00000205413__Unsupported
- Represents and autosomal dominant de novo

#### de_novo.vcf.bgz

- Similar to 1_labelled_variant (samples, call)
- Used in de novo classification test, so no Category flags

#### de_novo_ped.fam

- Pedigree for PROBAND, MOTHER, FATHER trio
- only PROBAND (male) is affected
- PLINK format

#### fake_gene_list.json

- 3 fake gene names as a JSON-encoded list
- 'foo bar', 'foo', 'bar'
- Used in testing the PanelApp methods

#### mock_pedigree.json

- JSON encoding members of 2 'families'
- Simulates format provided by SampleMetadata API
- Families 'FAM1' & 'FAM2'

#### mock_sm_lookup.json

- Mocked JSON lookup, mapping the sample IDs in mock_pedigree.json to a "CPGXXX" value
- List of Lists (representing SM-API's List of Tuples)

#### multiple_hail.vcf.bgz

- 2-sample VCF, containing a number of calls, annotated against different genes
- Used in the Hail Compound-Het calculation tests
- samples SAMPLE1, SAMPLE2

#### panel*.json

- collection of JSON values
- panelapp_* represents the results of PanelApp API calls
- panel_changes_expected - resultant values when two panels are compared
- panel_green... - result when panel filtered for 'green' genes only

#### pedfile.ped

- PED format file, contains 2 trios
- Used for a performance test of Peddy.Ped

#### seqr_tags.tsv

- Simulated data, representing the Seqr export format
- Contains 2 tagged variants
    - one has multiple tags, one is Talos-relevant
    - one has a single tag, not relevant

#### single_hail.vcf.bgz

- Single variant, single sample
- In the Hail classification tests, this is loaded as a MT
- Using pytest.parametrize, the same variant is loaded with different annotations, and run through the filtering process
- Sample ID = 'SAMPLE'

#### trio_plus_sibling.fam

- Plink format file, containing the sample family as de_novo_ped.fam
- Adds a second child, unaffected sibling

#### newphase.vcf.bgz  *(generated — not committed)*

- Two variants on chr20 at the same loci as 1_labelled_variant
- Both annotated against ENSG00000075043
- mother_1 carries both variants as phased hets (0|1) in phase-set PS=420
- male and father_1 are WT (0/0)
- Used to verify that same-phase-set compound-het pairs are correctly excluded
