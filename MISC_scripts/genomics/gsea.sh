java -cp /home1/zl222/software/gsea2.jar  -Xmx512m xtools.gsea.Gsea -res expression.txt -cls Pha4.cls#Log2Ratio_versus_Background -gmx Pha4.grp -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric Diff_of_Classes -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out output -gui false


