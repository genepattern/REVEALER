# REVEALER GpUnit tests
For now there is only one test, which corresponds exactly to simply running the module with its defaults.

While the test has no built-in file comparison, it's possible to compare the result PDF file at the 
command-line using ImageMagick:

    $ for (( i = 0; i < 6; i++ )); do \ 
        compare -metric RMSE BCAT_vs_MUT_CNA.pdf[$i] BETA_CATENIN_vs_MUT_CNA.pdf[$i] diff${i}.pdf; echo ""; done
    36.8224 (0.000561874)
    0 (0)
    0 (0)
    0 (0)
    0 (0)
    0 (0)

These ImageMagick calls are performing a page-by-page comparison (six pages in total) and will generate a 
difference-highlighting PDF for each page, diff0.pdf to diff5.pdf.  In the example comparison shown above, the
first page differs slightly (0.05% difference, as indicated by the parenthetical number).  The remaining pages
are identical.

In this example, the first page of the GpUnit output had some text in boldface that the expected result has 
in normal weight.  Checking the diff0.pdf difference-highlight shows no visually detectable differences, however,
and that the text itself is identical.  