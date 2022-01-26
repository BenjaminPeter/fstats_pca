rule svg2png:
    input: '{name}.svg'
    output: '{name}.png'
    shell: 'inkscape {input} --export-type=png --export-background="#fff" --export-dpi 600'
rule svg2pdf:
    input: '{name}.svg'
    output: '{name}.pdf'
    shell: 'inkscape {input} --export-type=pdf --export-background="#fff"'

rule allpng:
    input:
        'figures/fig_data_world.png',
        'figures/fig_data_europe.png',
        'figures/outgroup_f3.png',
        'figures/fstats_admixture_pca.png',
        'figures/fstats_pca_vs_tree.png'

rule allpdf:
    input:
        'figures/fig_data_world.pdf',
        'figures/fig_data_europe.pdf',
        'figures/outgroup_f3.pdf',
        'figures/fstats_admixture_pca.pdf',
        'figures/fstats_pca_vs_tree.pdf',
        'figures/fig_f4_ratio.pdf',
        'figures/fig_proj.pdf',

rule tables:
    input:
        'data/subdata/worldfoci2.ind',
        'data/subdata/westeurasian1.ind',
    output:
        'supp_tables/file1.txt',
        'supp_tables/file2.txt'
    shell:
        'cp {input[0]} {output[0]};'
        'cp {input[1]} {output[1]}'

rule all:
    input:
        rules.allpdf.input,
        rules.tables.output,
        tex='fstats_pca_lewontin.tex',
        bib='main.bib',
    output:
        'fstats_pca_lewontin.pdf'
    shell:
        'latexmk -pdf {input.tex}'

