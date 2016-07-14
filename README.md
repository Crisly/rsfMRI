# rsfMRI pipeline
This Nipype (Gorgolewski et al., 2011) pipeline does the preprocessing of rsFMRI data and consists of three main parts (functions):
* create_realign_flow
* create_reg_flow
* create_resting_preproc
    
The first bit is the realignment of 4-d data: "Create_realign_flow" realigns a time series to the middle volume using spline     interpolation and uses FSL MCFLIRT to realign the time series and ApplyWarp to apply the rigid
body transformations using spline interpolation.
    
The 2nd function registers the time series to the Montreal Neurological Institute (MNI) stereotaxic space by rigid-body          transformation using ANTs (Avants et al, 2008, Klein et al., 2009).
    
Finally a rsfMRI time series preprocessing workflow is used to filter out physiological noise using the motion-parameters obtained during realignment and the CSF-time course (CompCor; Behzadi, 2007). The data was resampled to a 3-mm isotropic resolution and moderately smoothed using an isotropic 5 mm at FWHM Gaussian kernel to improve signal to noise ratio.
    
#Links:
- Excellent tutorial: http://miykael.github.io/nipype-beginner-s-guide/index.html
- Main processing line: http://nipy.org/nipype/
- FSL tools: http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MCFLIRT
- Advances Normalization Tools (ANTs): http://stnava.github.io/ANTs/
- Paper from Klein et al., "ART, SyN (part of ANTs), IRTK, and SPM's DARTEL Toolbox gave the best results according to              overlap and distance measures, with ART and SyN delivering the most consistently high accuracy across subjects and label sets"; http://www.ncbi.nlm.nih.gov/pubmed/19195496
- **Lanting, C.**, Woźniak, A., van Dijk, P., & Langers, D. R. M. (2016). <a title="pdf" href="https://www.crisly.nl/docs/Lanting et al. - 2016 - Tinnitus- and Task-Related Differences in Resting-State Networks.pdf">Tinnitus- and Task-Related Differences in Resting-State Networks.</a> Advances in experimental medicine and biology, 894, 175-87. <a href="http://dx.doi.org/10.1007/978-3-319-25474-6_19">10.1007/978-3-319-25474-6_19</a>

#Graphs    
The various functione result in a graph describing the flow of the data:<br>
![realignment](https://github.com/Crisly/rsfMRI/blob/master/realign.png)
![normalization](https://github.com/Crisly/rsfMRI/blob/master/realign.png)

The complete graph looks like this:
![complete](https://github.com/Crisly/rsfMRI/blob/master/rsfMRI_complete.png)
    
