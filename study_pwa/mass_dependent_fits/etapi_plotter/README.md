Example command to draw wave projections, expects waves to be of the form {Lme}. For instance S_{0}^{+} is equal to S0++ and D2-+ is D_{-2}^{+}. L is flexible and can accomdate extra characters, for instance pD2-- can be used to represent a separate amplitude with the same angular dependence as D2-- but perhaps different mass dependence [i.e. a2(1320) and a2(1700)]

Flags
- `-s` semicolon separated string containng underscore separated strings to draw waves for. For instance S0++;D2++_pD2++ will request the S0++ wave to be drawn in one root file and another root file containing the interference between two waves: D2++ and pD2++
- `-a` Whether to acceptance correc the final output yields into the outputted log file
- `-var` Whether to plot all the variables or not, including angular variables for instance. Can save time if only interested in Mass plots
- `-gen` Whether to plot genmc weighted by amplitudes or not
- `-F` Amplitudes to output fit fractions for, outputted into the log file. Underscore separated. For instance S will output the fit fraction containing all the amplitudes that start with S.

Example Usage:

```
etapi_plotter etapi_result.fit -s S0+-;D2++_pD2++ -a true -var true -gen false -F S_D_pD
```


