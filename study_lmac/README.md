It turns out that sideband subtraction performs poorly near the etapi threshold in data. It seems to stem from omega to 3g and pi0pi0 to 4g combinatorial backgrounds.
Performs poorly = significant cusp like features appear where oversubtraction and undersubtraction seems to occur

The code in study_expectedYields is used for this study to generate apply the standard event selections, except the lmac cut of course. The state of the DSelector is saved in this folder.

To handle these backgrounds a new cut is introduced, denote low mass alternative combinations (lmac) removal. Typically a pi0 is paired from g1g2 and eta is formed from g3g4. Alternative combinations exist g1g3, g1g4, g2g3, g2g4. A box cut is made on low mass M(g1g3) vs M(g2g4) and M(g1g3) vs M(g2g3) which can be seen in the notebook.
