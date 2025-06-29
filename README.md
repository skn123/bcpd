
# Bayesian Coherent Point Drift / Domain Elastic Transform

This repository provides a software suite implementing Domain Elastic Transform (DET),
Bayesian Coherent Point Drift (BCPD/GBCPD), and Dependent Landmark Drift (DLD),
BCPD registers two point clouds, which can be applied to shape analysis, 3D model
reconstruction, and so forth. DET registers two functions, which can be applied
to aligning digital images and audio signals. DLD is a method for active shape model fitting.
All methods can be accelerated using downsampling and displacement field interpolation,
indicated by '++' after the method’s name, as in BCPD++.

For more information, see
[Hirose2022](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9918058) (GBCPD/GBCPD++),
[Hirose2020a](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8985307) (BCPD), and
[Hirose2020b](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9290402) (BCPD++).
Also, several examples can be watched in
[[Video 1]](https://youtu.be/OT97b60iBmQ),
[[Video 2]](https://youtu.be/pbLVMDj1Zro),
[[Video 3]](https://youtu.be/cET6gKAvjw0),
[[Video 4]](https://youtu.be/SoUTbH2tJj8).
If you have any questions, kindly email ohirose.univ+bcpd(at)gmail.com with your name and affiliation.

![alt text](https://github.com/ohirose/bcpd/blob/master/img/transfer.jpg?raw=true)

## Table of Contents

1. [References](#references)
2. [Performance](#performance)
    + [GBCPD vs CPD](#gbcpd-vs-cpd)
    + [BCPD vs CPD](#bcpd-vs-cpd)
    + [BCPD vs BCPD++](#bcpd-vs-bcpd)
3. [Demo](#demo)
    + [Dataset preparation](#dataset-preparation)
    + [Demo script execution](#demo-script-execution)
4. [Compilation](#compilation)
    + [Windows](#windows)
    + [MacOS and Linux](#macos-and-linux)
5. [Usage](#usage)
    + [Terms and symbols](#terms-and-symbols)
    + [Algorithms](#algorithms)
    + [Input data](#input-data)
    + [Tuning parameters](#tuning-parameters)
    + [Kernel functions](#kernel-functions)
6. [Acceleration](#acceleration)
    + [Nystrom method](#nystrom-method)
    + [KD tree search](#kd-tree-search)
    + [Downsampling](#downsampling)
    + [Interpolation](#interpolation)
6. [Options](#options)
    + [Convergence](#convergence)
    + [Normalization](#normalization)
    + [File output](#file-output)
    + [Terminal output](#terminal-output)

## References

The details of the algorithms are available in the following papers:
- [GBCPD/GBCPD++] O. Hirose,
  "[Geodesic-Based Bayesian Coherent Point Drift](https://ieeexplore.ieee.org/document/9918058),"
  [IEEE TPAMI](https://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=34), Oct 2022.
- [BCPD++] O. Hirose,
  "[Acceleration of non-rigid point set registration with downsampling and Gaussian process regression](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9290402),"
  [IEEE TPAMI](https://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=34), Dec 2020.
- [BCPD] O. Hirose,
  "[A Bayesian formulation of coherent point drift](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8985307),"
  [IEEE TPAMI](https://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=34), Feb 2020.
  - The article's [supplementary document](https://www.dropbox.com/s/pkgw2xxd0f3anfk/bcpd-appendix.pdf?dl=1)
    contains proofs of propositions.
  - Fig. 15 in the print version was accidentally replaced by Fig. 3 during the publication process after the review process. See an [erratum](https://ieeexplore.ieee.org/document/9506964) correcting the error.
  - In Fig. 2, following Proposition 3, please replace ![trace-before](https://github.com/ohirose/bcpd/blob/master/img/trace-before.png?raw=true) with ![trace-after](https://github.com/ohirose/bcpd/blob/master/img/trace-after.png?raw=true) as follows:

<img src="https://github.com/ohirose/bcpd/blob/master/img/correction.png" alt="correction" width="600"/>

## Performance

### GBCPD vs CPD
GBCPD works better than CPD and BCPD if it registers the shapes whose different parts neighboring each other:
![alt text](https://github.com/ohirose/bcpd/blob/master/img/gbcpd-cmp.png?raw=true)

### BCPD vs CPD
BCPD is faster than coherent point drift (CPD) and is often more accurate. The following figure shows a comparison
using Armadillo data included in demo data (vs Dr. Myronenko's implementation on Macbook Pro Early 2013):

<img src="https://github.com/ohirose/bcpd/blob/master/img/vs-cpd.png" alt="vs-cpd" width="400"/>

### BCPD vs BCPD++
BCPD++ is much faster than BCPD but is slightly less accurate (Mac Mini 2018).

<img src="https://github.com/ohirose/bcpd/blob/master/img/vs-plusplus.png" alt="vs-plusplus" width="400"/>

## Demo

The software provides the following demonstrations:
 - Function registration: image/shape/audio registration (DET)
 - Point set registration (BCPD/BCPD++)
 - 3D model reconstruction from point clouds (BCPD/BCPD++)
 - Shape transfer (BCPD/BCPD++)
 - Surface registration (GBCPD/GBCPD++)
 - Active shape model fitting (DLD)

### Dataset Preparation

Download the following datasets and move them into the `data` folder in this software:
- [3D reconstruction data](https://www.dropbox.com/scl/fi/j4ar6wmgixj7eare9yx5e/bcpd-data-3drecov-06Jul2025.zip?rlkey=m1z3r92uibv7d32dzyxyezcu7&dl=1).
- [DET data](https://www.dropbox.com/scl/fi/pqndryktjwzjxqmq9tcz5/det-demodata-10Jul2025.zip?rlkey=1kinhx9w06qwg21mulyv6cmc3&dl=1).
- [DLD data](https://www.dropbox.com/scl/fi/e9fbrxdkos9qej0enwkv1/dld-demodata06Jul2025.zip?rlkey=nu8snqlttc9fdmao8aql8l6r9&dl=1).
- [GBCPD data](https://www.dropbox.com/s/yssce2kmdil3fqs/gbcpd-demodata20220829.zip?dl=1).
- [BCPD++ data](https://www.dropbox.com/s/um46xujczko39jk/bcpd-pp-demodata20210226.zip?dl=1).
- [BCPD data](https://www.dropbox.com/s/6kd4uiyt150uyz9/bcpd-demodata20200127.zip?dl=1)

### Demo Script Execution:

All demo scripts except for the shape transfer can be executed as follows:

1. Start MATLAB and move into a foloder including demo scripts, e.g., `demo/bcpd-3drecov`.
2. Double-click a demo script, e.g., `demoRecovChef.m`.
3. Press the run button in the code editor of MATLAB.

The demo scripts of the shape transfer can be executed as follows:

1. Go to the `demo/shapeTransfer` folder using your terminal window.
2. Run a bash script, e.g., type `./shapeTransferA.sh` in the terminal.
3. Check output files named `transferV[1/2]_y.interpolated.obj`.

## Compilation

Please contact the author if the compilation fails.

### Windows

Ready to go. The compilation is not required. Use the binary file `bcpd.exe` in the `win` directory.

### macOS

1. Install Xcode, Xcode command line tools, Homebrew.
2. Install OpenMP and OpenBLAS using Homebrew.
3. Download and decompress the zip file containing bcpd source codes.
4. Move into the top directory of the uncompressed folder using the terminal window.
5. Type `make`.

### Linux

1. Install OpenMP and OpenBLAS.
2. Download and decompress the zip file containing bcpd source codes.
3. Move into the top directory of the uncompressed folder using the terminal window.
4. Type `make`.

## Usage

Brief instructions are printed by typing `./bcpd -v` (or `bcpd -v` for windows) in the terminal window.
The binary file can also be executed using the `system` function in MATLAB.
See MATLAB scripts in the `demo` folder regarding the usage of the binary file.

### Terms and symbols

- X: Target point set. The point set corresponding to the reference shape.
- Y: Source point set. The point set to be deformed. The mth point in Y is denoted by ym.
- N: The number of points in the target point set.
- M: The number of points in the source point set.
- D: Dimensionality of the space in which the source and target point sets are embedded.
- D': Dimensionality of the codomain used only for function registration.
- Fx: The matrix collecting the target function values.
- Fy: The matrix collecting the source function values.

### Algorithms

**[BCPD]** Type the following command in the terminal window for Mac/Linux:

- ` ./bcpd -x <target:X> -y <source:Y> (+options) `

**[DET]** Type the following command in the terminal window for Mac/Linux:

- ` ./bcpd -x <target:X> -X <target:Fx> -y <source:Y> -Y <source:Fy> (+options) `

**[DLD]** Type the following command in the terminal window for Mac/Linux:

- ` ./bcpd -x <target:X>  -y <mean:Y> -C <shape variations> (+options) `

For Windows, replace `./bcpd` with `bcpd` and run the command in the DOS prompt.

### Input data

- `-x [file]`: The target shape represented as a matrix of size N x D.
- `-y [file]`: The source shape represented as a matrix of size M x D.
- `-X [file]`: The target function values represented as Fx of size N x D' (DET only).
- `-Y [file]`: The source function values represented as Fy of size M x D' (DET only).
- `-C [file]`: The shape covariance's eigenvalues `L` and eigenvectors `Q`, formatted as `[L';Q]` (DLD only).

If the file names of target and source point sets are `X.txt` and `Y.txt`, these arguments can be omitted.

### Transformation models

- `-Tsrn`: Default. Similarity + nonrigid transformation, defined as T(y)=sR(y+v)+t.
- `-Tan`: Affine + nonrigid transformation, defined as T(y)=A(y+v)+t.
- `-Ta`: Affine transformation, defined as T(y)=Ay+t.
- `-Tsr`: Similarity transformation, defined as T(y)=sRy+t.
- `-Tr`: Rigid transformation, defined as T(y)=Ry+t.
- `-Tn`: Nonrigid transformation, defined as T(y)=y+v.

Note s is a scale factor, R is a rotation matrix, A is a matrix of full rank,
t is a translation vector, and v is a pointwise displacement.

### Tuning parameters

- `-w [real]`: Omega. Outlier probability in (0,1).
- `-l [real]`: Lambda. Positive. It controls the expected length of deformation vectors. Smaller is longer.
- `-b [real]`: Beta. Positive. It controls the range where deformation vectors are smoothed.
- `-g [real]`: Gamma. Positive. It controls how much the initial alignment are ignored.
- `-k [real]`: Kappa. Positive. It controls the randomness of mixing coefficients.

The expected length of deformation vectors is sqrt(D/lambda). Set gamma around 2 to 10 if your target point set
is largely rotated. If input shapes are roughly registered, use `-g0.1` with the option `-ux`.
The default kappa is infinity, which means that all mixing coefficients are equivalent.
Do not specify `-k infinity` or extremely large kappa to impose the equivalence of mixing coefficients,
which sometimes causes an error.

### Kernel functions

#### Standard kernels:

- `-G [1-3]`: Switch kernel functions.
  - `-G1` Inverse multiquadric: `(||ym-ym'||^2+beta^2)^(-1/2)`
  - `-G2` Rational quadratic: `1-||ym-ym'||^2/(||ym-ym'||^2+beta^2)`
  - `-G3` Laplace: `exp(-|ym-ym'|/beta)`

The Gaussian kernel `exp(-||ym-ym'||^2/2*beta^2)` is used unless specified.
Here, `ym` represents the mth point in Y. The tuning parameter of a kernel functions is denoted by beta,
which controls the range where deformation vectors are smoothed.

#### Geodesic kernel:

- `-G [string,real,file]`: Geodesic kernel with an input mesh. E.g., `-G geodesic,0.2,triangles.txt`.
  - 1st argument: The string `geodesic` only., i.e., the tag representing the geodesic kernel.
  - 2nd argument: Tau. The rate controlling the balance between geodesic and Gaussian kernels.
  - 3rd argument: The file that defines a triangle mesh.

- `-G [string,real,int,real]`: Geodesic kernel without an input mesh. E.g., `-G geodesic,0.2,8,0.15`.
  - 1st argument: The string `geodesic` only, i.e., the tag representing the geodesic kernel.
  - 2nd argument: Tau. The rate controlling the balance between geodesic and Gaussian kernels.
  - 3rd argument: The number of neighbors for each node, required for k-NN graph construction.
  - 4th argument: The radius that defines neighbors for each node, required for k-NN graph construction.

The geodesic kernel usually outperforms standard kernels when different parts of a source shape are closely located.
If the mesh of a source shape is available, choose the first option; it typically works better than the second option.
Otherwise, choose the second option; BCPD automatically creates the graph required for geodesic computations.
The mesh file must be a tab-separated file that contains three integers for each line; a triangle is defined
as a triplet of vertices. The following parameters tune the geodesic kernel:

- `-b [real]`: Beta. Positive. Gaussian function's width.
- `-K [int]`:  K. Positive. Rank constraint on G.
- `-z [real]`: Epsilon. Positive. Acceptable condition number of G.

## Acceleration
![alt text](https://github.com/ohirose/bcpd/blob/master/img/lucy.png?raw=true)

BCPD/DLD/DET can be accelerated inside and outside variational Bayes inference (VB), separately. The Nystrom method and
KD-tree search accelerate VBI. The former works before approaching convergence, whereas the latter works near
convergence. The following option accelerates VBI with default parameters:

- `-A`: VB acceleration with parameters, i.e., `-K70 -J300 -p -d7 -e0.15 -f0.2`.

Downsampling and deformation vector interpolation, called BCPD++/DET++/DLD++, accelerate non-rigid registration
outside VBI. For example, the following options activate the acceleration:

- `-DB,5000,0.08`: Downsampling/Interporation acceleration outside VBI.

If N and M are larger than several thousand, activate either the internal or external acceleration.
If N and M are more than several hundreds of thousands, activate both accelerations.
Also, the acceleration methods reduce memory consumption. Either internal or external acceleration method
**MUST** be activated to avoid a memory allocation error for such a large dataset.
Otherwise, the software sometimes fails without any error notice.

### Nystrom method

- `-K [int]`: #Nystrom samples for computing G.
- `-J [int]`: #Nystrom samples for computing P.
- `-r [int]`: Random number seed for the Nystrom method. Reproducibility is guaranteed if the same number is specified.

Specify `-J300 -K70`, for example. The acceleration based only on the Nystrom method probably fails to converge;
do not forget activating the KD-tree search. Disable `-J` when performing function registration to avoid unstable computations.

### KD tree search

- `-p`: KD-tree search is turned on if specified. The following options fine-tune the KD tree search.
  - `-d [real]`: Scale factor of sigma that defines areas to search for neighbors.
  - `-e [real]`: Maximum radius to search for neighbors.
  - `-f [real]`: The value of sigma at which the KD tree search is turned on.

The default parameters are `-d7 -e0.15 -f0.2`.
Retry the execution with `-p -f0.3` unless the Nystrom method is replaced by the KD-tree search during optimization.

### Downsampling

- `-D [char,int,real]`: Changes the number of points. E.g., `-D'B,10000,0.03'`.
  - 1st argument: One of the symbols: [X,Y,B,x,y,b]; x: target; y: source; b: both, upper: voxel, lower: ball.
  - 2nd argument: The number of points to be extracted by the downsampling.
  - 3rd argument: The voxel size or ball radius required for downsampling.

Input point sets can be downsampled by the following algorithms:
1. voxel-grid resampling with voxel width r,
2. ball resampling with the radius r, and
3. random resampling with equivalent sampling probabilities.

The parameter r can be specified as the 3rd argument of `-D`. If r is specified as 0,
sampling scheme (3) is selected. The numbers of points in downsampled target and source point sets
can be different; specify the `-D` option twice, e.g., `-D'X,6000,0.08' -D'Y,5000,0.05'`.
For more information, see [paper](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9290402) and
[appendix](https://ieeexplore.ieee.org/ielx7/34/4359286/9290402/supp1-3043769.pdf?tp=&arnumber=9290402).

### Interpolation

Downsampling automatically activates the deformation vector interpolation.
The resulting registered shape with interpolation is output to the file with the suffix `y.txt`.

## Options

Default values will be used unless specified.

### Convergence

- `-c [real]`: Convergence tolerance.
- `-n [int ]`: The maximum number of VB loops.
- `-N [int ]`: The minimum number of VB loops.

The default minimum VB iteration is `30`, which sometimes causes an error for small data.
If the bcpd execution stopped within 30 loops with an error notice, execute it again after
setting `-N1`, which removes the constraint on the minimum VB iteration.
The default value of the convergence tolerance is `1e-4`. If your point sets are smooth
surfaces with moderate numbers of points, specify `-c 1e-5` or `-c 1e-6`.

### Normalization

- `-u [char]`: Chooses a normalization option by specifying the argument of the option, e.g., `-ux`.
  - `e`: Each of X and Y is normalized separately (default).
  - `x`: X and Y are normalized using the location and the scale of X.
  - `y`: X and Y are normalized using the location and the scale of Y.
  - `n` : Normalization is skipped (not recommended).

Using `-ux` or `-uy` is recommended with `-g0.1` if input point sets are roughly registered.
The option `-un` is not recommended because choosing beta and lambda becomes non-intuitive.

### Terminal output

- `-v`: Print the version and the simple instruction of this software.
- `-q`: Quiet mode. Print nothing.
- `-W`: Disable warnings.
- `-h`: History mode. Alternative terminal output regarding optimization.

### File output

- `-o [string]`: Prefix of output file names.
- `-s [string]`: Save variables by specifying them as the argument of the option, e.g., `-sYP`.
  - `y`: Resulting deformed shape (=y).
  - `x`: Target shape with alignment (=x).
  - `v`: Displacement vector (=v).
  - `c`: non-outlier labels (=c).
  - `e`: matched points (=e).
  - `a`: Mixing coefficients (=alpha).
  - `P`: Nonzero matching probabilities (=P).
  - `T`: Similarity transformation (=s,R,t).
  - `Y`: Optimization trajectory.
  - `t`: Computing time (real/cpu) and sigma for each loop.
  - `A`: All of the above.

The resulting deformed shape y will be output without the `-s` option.
If `Y` is specified as an argument of `-s`, the optimization trajectory
will be saved to the binary file `.optpath.bin`.
The trajectory can be viewed using scripts: `optpath.m` for 2D data and
`optpath3.m` for 3D data. Saving a trajectory is memory-inefficient. Disable it if both N and M
are more than several hundreds of thousands. If `P` is specified as an argument of `-s`,
nonzero elements of matching probability P will be output.

