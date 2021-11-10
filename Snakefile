rule svg2png:
    input: '{name}.svg'
    output: '{name}.png'
    shell: 'inkscape {input} --export-type=png --export-background="#fff" --export-dpi 600'
rule svg2pdf:
    input: '{name}.svg'
    output: '{name}.pdf'
    shell: 'inkscape {input} --export-type=pdf --export-background="#fff"'

rule all:
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

rule latex:
    input:
        rules.allpdf.input,
        tex='fstats_pca_lewontin.tex'
    output:
        'fstats_pca_lewontin.pdf'
    shell:
        'latexmk -pdf {input.tex}'
