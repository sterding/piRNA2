for i in `ls -d ../mouse_*dpp_wt_smallRNA_ox_rep1`;
do
    echo $i;
    samplename=${i/..\//}
    cd $i
    sum__readsall_speciesall_readsuniq_speciesuniq=`awk '{split($4,a,"_");if($5==1) uniq[$4]=a[2]; all[$4]=a[2];}END{n=0;s=0;n2=0;s2=0;for(i in all) {n++;s=s+all[i]}; for(i in uniq) {n2++;s2=s2+uniq[i];} print s,n, s2, n2;}' allmap.bed_2_piRNA`
    Rscript $HOME/projects/piRNA/src/plot_quantify_smallRNA genes.allmap.nointron.ppm.rpkm genes.uniqmap.nointron.ppm.rpkm $samplename.gene.nointron.piRNA $sum__readsall_speciesall_readsuniq_speciesuniq
    cd ..
done
