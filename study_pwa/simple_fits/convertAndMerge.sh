for t in 010020 0200325 0325050 050075 075100
do
    convert -define pdf:use-cropbox=true -density 300 -quality 100 -sharpen 0x1.0 results/bw_plus_poly_t$t.pdf results/bw_plus_poly_t$t.png
    convert -define pdf:use-cropbox=true -density 300 -quality 100 -sharpen 0x1.0 results/voigt_plus_poly_t$t.pdf results/voigt_plus_poly_t$t.png
done

montage results/bw*png -geometry +8+8 -tile 2x3 -quality 90 results/montage_bw.png
montage results/voigt*png -geometry +8+8 -tile 2x3 -quality 90 results/montage_voigt.png
