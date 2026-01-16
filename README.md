Compile:
```
gfortran Index.f90  MODULES.f90 -o ./a.out
```
Run:
```
./a.out
```
Plot the UV map
```
python plot_uv_map.py --input Output/E_Calc_12Dec.UV --grid Input/grid --output uvi_map_interpolated.png
```