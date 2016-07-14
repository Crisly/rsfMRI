import os                                    # system functions
import os.path as op
import numpy as np
import nipype.interfaces.io as nio           # Data i/o
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl          # fsl
from nipype.algorithms.misc import TSNR
import nibabel as nb
from nipype.interfaces.utility import Function
fsl.FSLCommand.set_default_output_type('NIFTI_GZ')
#from nipype import config
#config.enable_debug_mode()
""" This Nipype (Gorgolewski et al., 2011) function does the preprocessing of rsFMRI data and consists of three main parts (functions):
    - create_realign_flow
    - create_reg_flow
    - create_resting_preproc
    
    First some important bits are defined about the project: which anatomy files do you want (FFE or TFE),
    which anatomy template, what TR, fwhm, filter-settings etc.
    
    The second bit (first main function) is the realignment of 4-d data: "Create_realign_flow" realigns a time series to the         middle volume using spline interpolation and uses FSL MCFLIRT to realign the time series and ApplyWarp to apply the rigid
    body transformations using spline interpolation.
    
    The third bit (2nd function) registers the time series to the Montreal Neurological Institute (MNI) stereotaxic space by         rigid-body transformation using ANTs (Avants et al, 2008, Klein et al., 2009).
    
    Finally a "resting" time series preprocessing workflow is used to filter out physiological noise using the motion-parameters     obtained during realignment and the CSF-time course (CompCor; Behzadi, 2007). The data was resampled to a 3-mm isotropic         resolution and moderately smoothed using an isotropic 5 mm at FWHM Gaussian kernel to improve signal to noise ratio.
    
    Links:
    - Excellent tutorial: http://miykael.github.io/nipype-beginner-s-guide/index.html
    - Main processing line: http://nipy.org/nipype/
    - FSL tools: http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MCFLIRT
    - Advances Normalization Tools (ANTs): http://stnava.github.io/ANTs/
    - Paper from Klein et al., "ART, SyN (part of ANTs), IRTK, and SPM's DARTEL Toolbox gave the best results according to              overlap and distance measures, with ART and SyN delivering the most consistently high accuracy across subjects and label         sets"; http://www.ncbi.nlm.nih.gov/pubmed/19195496
    - finally: see Lanting, C., Wozniak, A., van Dijk, P., & Langers, D. R. M. (2016). "Tinnitus- and Task-Related                      Differences in Resting-State Networks." Advances in experimental medicine and biology, 894, 175-87. doi:                         10.1007/978-3-319-25474-6_19; http://link.springer.com/chapter/10.1007%2F978-3-319-25474-6_19
""" 

def remove_interleaved(in_file):
    """Function to remove the image artefacts (along the z-direciton) created if your data is acquired in an interleaved fashion
    Parameters
    ----------
    in_files: one Nifti 4D time series
    
    Returns
    -------
    out_file: one Nifti 4D time series
    """
    import os
    import nibabel as nb
    print 'read: ' + str(in_file)
    path, filename = os.path.split(in_file)
    img = nb.load(in_file)
    data = img.get_data()
    affine = img.get_affine()
    deinterlaced = (data[:,:,0:-1,:] + data[:,:,1:,:])/2
    deinterlaced_file = nb.Nifti1Image(deinterlaced, affine)
    newfile = os.path.join(path,'z'+filename)
    nb.save(deinterlaced_file, newfile)
    print 'wrote: ' + newfile
    return newfile,deinterlaced

def median(in_files):
    """Computes an average of the median of each realigned timeseries
    Parameters
    ----------
    in_files: one or more realigned Nifti 4D time series
    
    Returns
    -------
    out_file: a 3D Nifti file
    """
    average = None
    for idx, filename in enumerate(filename_to_list(in_files)):
        img = nb.load(filename)
        data = np.median(img.get_data(), axis=3)
        if not average:
            average = data
        else:
            average = average + data
    median_img = nb.Nifti1Image(average/float(idx + 1),
                                img.get_affine(), img.get_header())
    filename = os.path.join(os.getcwd(), 'median.nii.gz')
    median_img.to_filename(filename)
    return filename

def get_substitutions(subject_id):
    '''Replace output names of files with more meaningful ones
    '''
    return [('vol0000_warp_merged_detrended_regfilt_filt_trans_mask_smooth',
             '%s_filtered_smooth'%subject_id),
            ('vol0000_warp_merged_tsnr_stddev_thresh',
             '%s_noisyvoxels'%subject_id)]

def select_volume(filename, which):
    """Return the middle index of a file
    """
    from nibabel import load
    import numpy as np
    if which.lower() == 'first':
        idx = 0
    elif which.lower() == 'middle':
        idx = int(np.ceil(load(filename).get_shape()[3]/2))
    else:
        raise Exception('unknown value for volume selection : %s'%which)
    return idx

def extract_noise_components(realigned_file, noise_mask_file, num_components):
    """ Derive components most reflective of physiological noise
    """
    import os
    from nibabel import load
    import numpy as np
    import scipy as sp
    from scipy.signal import detrend
    imgseries = load(realigned_file)
    noise_mask = load(noise_mask_file)
    voxel_timecourses = imgseries.get_data()[np.nonzero(noise_mask.get_data())]
    for timecourse in voxel_timecourses:
        timecourse[:] = detrend(timecourse, type='constant')
    u,s,v = sp.linalg.svd(voxel_timecourses, full_matrices=False)
    components_file = os.path.join(os.getcwd(), 'noise_components.txt')
    np.savetxt(components_file, v[:num_components, :].T)
    return components_file

def create_realign_flow(name='realign'):
    """Realign a time series to the middle volume using spline interpolation
    Uses MCFLIRT to realign the time series and ApplyWarp to apply the rigid
    body transformations using spline interpolation (unknown order).
    
    Example
    -------
    >>> wf = create_realign_flow()
    >>> wf.inputs.inputspec.func = 'f3.nii'
    >>> wf.run() # doctest: +SKIP
    """
    realignflow = pe.Workflow(name=name)
    
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['func',
                                                                 ]),
                        name='inputspec')
    outputnode = pe.Node(interface=util.IdentityInterface(fields=[
                                                               'realigned_file',
                                                                 ]),
                        name='outputspec')
    realigner = pe.Node(fsl.MCFLIRT(save_mats=True, stats_imgs=True),
                        name='realigner')
    splitter = pe.Node(fsl.Split(dimension='t'), name='splitter')
    warper = pe.MapNode(fsl.ApplyWarp(interp='spline'),
                        iterfield=['in_file', 'premat'],
                        name='warper')
    joiner = pe.Node(fsl.Merge(dimension='t'), name='joiner')
    
    # remove zebras before realignment in case of RS-data
    deinterlacer = pe.Node(interface=Function(input_names=['in_file', ],
                        output_names=['deinterlaced','data' ],
                        function=remove_interleaved), name='deinterlacer')
    realignflow.connect(inputnode, 'func', deinterlacer,'in_file')
    realignflow.connect(deinterlacer, 'deinterlaced', realigner, 'in_file')
    realignflow.connect(deinterlacer, ('deinterlaced', select_volume, 'middle'),
                        realigner, 'ref_vol')
    realignflow.connect(realigner, 'out_file', splitter, 'in_file')
    realignflow.connect(realigner, 'mat_file', warper, 'premat')
    realignflow.connect(realigner, 'variance_img', warper, 'ref_file')
    realignflow.connect(splitter, 'out_files', warper, 'in_file')
    realignflow.connect(warper, 'out_file', joiner, 'in_files')
    realignflow.connect(joiner, 'merged_file', outputnode, 'realigned_file')
    realignflow.write_graph(dotfilename="realign",graph2use='colored', format='png', simple_form=True)
    return realignflow


def create_reg_flow(name='registration'):
    from nipype.interfaces.ants import Registration
    from nipype.interfaces.ants import ApplyTransforms
    from nipype.interfaces.utility import Merge
    import nipype.pipeline.engine as pe
    
    """Register time series to template through subjects' T1
    Example
    -------
    >>> wf.inputs.inputspec.func = op.join(resdir,'func.nii.gz')
    >>> wf.inputs.inputspec.struct = op.join(resdir,'t1.nii.gz')
    >>> wf.inputs.inputspec.Template = op.join(resdir,'MNI152_T1_2mm_brain.nii.gz')
    >>> wf.inputs.inputspec.Template_3mm = op.join(resdir,'MNI152_T1_3mm_brain.nii.gz')
    >>> wf.run()
    """
    def select_volume(filename, which):
        """Return the middle index of a file
        """
        from nibabel import load
        import numpy as np
        if which.lower() == 'first':
            idx = 0
        elif which.lower() == 'middle':
            idx = int(np.ceil(load(filename).get_shape()[3]/2))
        else:
            raise Exception('unknown value for volume selection : %s'%which)
        return idx
    
    extract_ref = pe.Node(interface=fsl.ExtractROI(t_size=1),
                      name = 'extractref')
    
    # first 'bet' the anatomy file, including bias reduction step
    bet = pe.Node(interface=fsl.BET(),mask=True, name='bet')
    bet.inputs.reduce_bias = True #or false?!
    
    # coregistration step based on rigid transformation using ANTs
    # see e.g., # see https://github.com/binarybottle/mindboggle/issues/15
    
    coreg = pe.Node(Registration(), name='CoregAnts')
    coreg.inputs.output_warped_image = 'func2highres_.nii.gz'
    coreg.inputs.output_transform_prefix = "func2highres_"
    coreg.inputs.transforms = ['Rigid']
    coreg.inputs.transform_parameters = [(0.1,), (0.1,)]
    coreg.inputs.number_of_iterations = [[1000,500,250,100]]
    coreg.inputs.dimension = 3
    coreg.inputs.write_composite_transform = True
    coreg.inputs.collapse_output_transforms = True
    coreg.inputs.metric = ['MI']
    coreg.inputs.metric_weight = [1]
    coreg.inputs.radius_or_number_of_bins = [32]
    coreg.inputs.sampling_strategy = ['Regular']
    coreg.inputs.sampling_percentage = [0.25]
    coreg.inputs.convergence_threshold = [1.e-8]
    coreg.inputs.convergence_window_size = [10]
    coreg.inputs.smoothing_sigmas = [[3,2,1,0]]
    coreg.inputs.sigma_units = ['mm']
    coreg.inputs.shrink_factors = [[8,4,2,1]]
    coreg.inputs.use_estimate_learning_rate_once = [True]
    coreg.inputs.use_histogram_matching = [False]
    coreg.inputs.initial_moving_transform_com = True
    coreg.inputs.output_warped_image = True
    coreg.inputs.winsorize_lower_quantile = 0.01
    coreg.inputs.winsorize_upper_quantile = 0.99
    
    # registration or normalization step based on symmetric diffeomorphic image registration (SyN) using ANTs
    # see for parameters: https://gist.github.com/satra/8439778
    
    reg = pe.Node(Registration(), name='NormalizationAnts')
    reg.inputs.output_transform_prefix = 'highres2template'
    reg.inputs.output_warped_image = 'highres2template.nii.gz'
    
    reg.inputs.output_transform_prefix = "highres2template_"
    reg.inputs.transforms = ['Rigid', 'Affine', 'SyN']
    reg.inputs.transform_parameters = [(0.1,), (0.1,), (0.2, 3.0, 0.0)]
    # for testing; means less number of iterations
    # reg.inputs.number_of_iterations = ([[10000, 111110, 11110]] * 2 + [[10, 5, 2]])
    
    # final; more iterations
    reg.inputs.number_of_iterations = ([[10000, 111110, 11110]] * 2 + [[100, 50, 30]])
    reg.inputs.dimension = 3
    reg.inputs.write_composite_transform = True
    reg.inputs.collapse_output_transforms = True
    reg.inputs.initial_moving_transform_com = True
    reg.inputs.metric = ['Mattes'] * 2 + [['Mattes', 'CC']]
    reg.inputs.metric_weight = [1] * 2 + [[0.5, 0.5]]
    reg.inputs.radius_or_number_of_bins = [32] * 2 + [[32, 4]]
    reg.inputs.sampling_strategy = ['Regular'] * 2 + [[None, None]]
    reg.inputs.sampling_percentage = [0.3] * 2 + [[None, None]]
    reg.inputs.convergence_threshold = [1.e-8] * 2 + [-0.01]
    reg.inputs.convergence_window_size = [20] * 2 + [5]
    reg.inputs.smoothing_sigmas = [[4, 2, 1]] * 2 + [[1, 0.5, 0]]
    reg.inputs.sigma_units = ['vox'] * 3
    reg.inputs.shrink_factors = [[3, 2, 1]]*2 + [[4, 2, 1]]
    reg.inputs.use_estimate_learning_rate_once = [True] * 3
    reg.inputs.use_histogram_matching = [False] * 2 + [True]
    reg.inputs.winsorize_lower_quantile = 0.005
    reg.inputs.winsorize_upper_quantile = 0.995
    reg.inputs.args = '--float'
    
    # fetch input
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['func',
                                                                 'struct',
                                                                 'Template',
                                                                 'Template_3mm'
                                                                 ]),
                        name='inputspec')
    
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['registered_func','registered_T1']),
                        name='outputspec')
    # combine transforms
    pickfirst = lambda x: x[0]
    merge = pe.MapNode(Merge(2), iterfield=['in2'], name='mergexfm')
    
    # apply the combined transform
    applyTransFunc = pe.MapNode(ApplyTransforms(),iterfield=['input_image', 'transforms'],
                         name='applyTransFunc')
    applyTransFunc.inputs.input_image_type = 3
    applyTransFunc.inputs.interpolation = 'BSpline'
    applyTransFunc.inputs.invert_transform_flags = [False, False]
    applyTransFunc.inputs.terminal_output = 'file'
    
    regworkflow = pe.Workflow(name=name)
    regworkflow.connect(inputnode, 'struct', bet, 'in_file')
    regworkflow.connect(bet,'out_file', coreg, 'fixed_image')
    regworkflow.connect(inputnode, 'func', extract_ref, 'in_file')
    regworkflow.connect(inputnode, ('func', select_volume, 'middle'), extract_ref, 't_min')
    regworkflow.connect(extract_ref, 'roi_file', coreg, 'moving_image')
    regworkflow.connect(bet, 'out_file', reg, 'moving_image')
    regworkflow.connect(inputnode, 'Template', reg, 'fixed_image')
    
    # get transform of functional image to template and apply it to the functional images
    # to template_3mm (same space as template)
    # regworkflow.connect(inputnode, 'Template_3mm', applyTransFunc, 'reference_image')
    # alternatively, just use the overall template
    
    regworkflow.connect(inputnode, 'Template', applyTransFunc, 'reference_image')
    regworkflow.connect(inputnode, 'func', applyTransFunc, 'input_image')
    regworkflow.connect(coreg, ('composite_transform', pickfirst), merge, 'in1')
    regworkflow.connect(reg, ('composite_transform', pickfirst), merge, 'in2')
    regworkflow.connect(merge, 'out', applyTransFunc, 'transforms')
    
    #output
    regworkflow.connect(applyTransFunc, 'output_image', outputnode, 'registered_func')
    regworkflow.connect(reg, 'warped_image',outputnode, 'registered_T1')
    regworkflow.write_graph(dotfilename="normalization",graph2use='colored', format='png', simple_form=True)
    return regworkflow

def create_resting_preproc(name='restpreproc'):
    """Create a "resting" time series preprocessing workflow
    The noise removal is based on Behzadi et al. (2007)
    Parameters
    ----------
    name : name of workflow (default: restpreproc)
    Inputs::
        inputspec.func : functional run (filename or list of filenames)
    Outputs::
        outputspec.noise_mask_file : voxels used for PCA to derive noise components
        outputspec.filtered_file : bandpass filtered and noise-reduced time series
    Example
    -------
    >>> TR = 3.0
    >>> wf = create_resting_preproc()
    >>> wf.inputs.inputspec.func = 'f3.nii'
    >>> wf.inputs.inputspec.num_noise_components = 6
    >>> wf.inputs.inputspec.highpass_sigma = 100/(2*TR)
    >>> wf.inputs.inputspec.lowpass_sigma = 12.5/(2*TR)
    >>> wf.run() # doctest: +SKIP
    """
    restpreproc = pe.Workflow(name=name)
    # Define nodes
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['func',
                                                                 'struct',
                                                                 'template',
                                                                 'template_3mm',
                                                                 'num_noise_components',
                                                                 'highpass_sigma',
                                                                 'lowpass_sigma',
                                                                 'fwhm'
                                                                 ]),
                        name='inputspec')
    outputnode = pe.Node(interface=util.IdentityInterface(fields=[
                                                              'noise_mask_file',
                                                              'filtered_file',
                                                              ]),
                     name='outputspec')
    pickfirst = lambda x: x[0]
    realigner = create_realign_flow()
    registration = create_reg_flow()
    
    tsnr = pe.Node(TSNR(regress_poly=2), name='tsnr')
    getthresh = pe.Node(interface=fsl.ImageStats(op_string='-p 98'),
                           name='getthreshold')
    threshold_stddev = pe.Node(fsl.Threshold(), name='threshold')
    compcor = pe.Node(util.Function(input_names=['realigned_file',
                                                 'noise_mask_file',
                                                 'num_components'],
                                     output_names=['noise_components'],
                                     function=extract_noise_components),
                       name='compcorr')
    remove_noise = pe.Node(fsl.FilterRegressor(filter_all=True),
                           name='remove_noise')
    bandpass_filter = pe.Node(fsl.TemporalFilter(),
                              name='bandpass_filter')
    #smooth data
    bet2 = pe.Node(interface=fsl.BET(),mask=True, name='bet2')
    medianval = pe.MapNode(interface=fsl.ImageStats(op_string='-k %s -p 50'),
                       iterfield = ['in_file'],
                       name='medianval')
    meanfunc2 = pe.MapNode(interface=fsl.ImageMaths(op_string='-Tmean',
                                                suffix='_mean'),
                       iterfield=['in_file'],
                       name='meanfunc2')
    maskfunc2 = pe.MapNode(interface=fsl.ImageMaths(suffix='_mask',
                                                op_string='-mas'),
                       iterfield=['in_file'],
                       name='maskfunc2')
    mergenode = pe.Node(interface=util.Merge(2, axis='hstack'),
                    name='merge')
    smooth = pe.MapNode(interface=fsl.SUSAN(),
                    iterfield=['in_file', 'brightness_threshold','usans',],
                    name='smooth')
    smooth.inputs.fwhm = fwhm
    
    # Define connections
    restpreproc.connect(inputnode, 'func', realigner, 'inputspec.func')
    restpreproc.connect(realigner, 'outputspec.realigned_file', tsnr, 'in_file')
    restpreproc.connect(tsnr, 'stddev_file', threshold_stddev, 'in_file')
    restpreproc.connect(tsnr, 'stddev_file', getthresh, 'in_file')
    restpreproc.connect(getthresh, 'out_stat', threshold_stddev, 'thresh')
    restpreproc.connect(realigner, 'outputspec.realigned_file',
                        compcor, 'realigned_file')
    restpreproc.connect(threshold_stddev, 'out_file',
                        compcor, 'noise_mask_file')
    restpreproc.connect(inputnode, 'num_noise_components',
                        compcor, 'num_components')
    restpreproc.connect(tsnr, 'detrended_file',
                        remove_noise, 'in_file')
    restpreproc.connect(compcor, 'noise_components',
                        remove_noise, 'design_file')
    restpreproc.connect(inputnode, 'highpass_sigma',
                        bandpass_filter, 'highpass_sigma')
    restpreproc.connect(inputnode, 'lowpass_sigma',
                        bandpass_filter, 'lowpass_sigma')
    restpreproc.connect(remove_noise, 'out_file', bandpass_filter, 'in_file')
    restpreproc.connect(threshold_stddev, 'out_file',
                        outputnode, 'noise_mask_file')
    restpreproc.connect(bandpass_filter, 'out_file',registration, 'inputspec.func')
    restpreproc.connect(inputnode, 'struct', registration, 'inputspec.struct')
    restpreproc.connect(inputnode, 'template', registration, 'inputspec.Template')
    restpreproc.connect(inputnode, 'template_3mm', registration, 'inputspec.Template_3mm')
    restpreproc.connect(registration, 'outputspec.registered_func', medianval, 'in_file')
    
    #select only first of output to bet2
    restpreproc.connect(registration, ('outputspec.registered_func',pickfirst), bet2, 'in_file')
    restpreproc.connect(bet2, 'out_file', medianval, 'mask_file')
    restpreproc.connect(registration, 'outputspec.registered_func', maskfunc2, 'in_file')
    restpreproc.connect(bet2, 'out_file', maskfunc2, 'in_file2')
    restpreproc.connect(maskfunc2, 'out_file', meanfunc2, 'in_file')
    restpreproc.connect(meanfunc2,'out_file', mergenode, 'in1')
    restpreproc.connect(medianval,'out_stat', mergenode, 'in2')
    
    def getbtthresh(medianvals):
        return [0.75*val for val in medianvals]
    
    def getusans(x):
        return [[tuple([val[0],0.75*val[1]])] for val in x]
    
    restpreproc.connect(maskfunc2, 'out_file', smooth, 'in_file')
    restpreproc.connect(medianval, ('out_stat', getbtthresh), smooth, 'brightness_threshold')
    restpreproc.connect(mergenode, ('out', getusans), smooth, 'usans')
    restpreproc.connect(smooth, 'smoothed_file', outputnode,'filtered_file')
    restpreproc.write_graph(dotfilename="rsfMRI_complete",graph2use='colored', format='png', simple_form=True)
    return restpreproc
    
    
## important information goes here ##
resdir = '/Users/cplanting/ProjectsScience/rsfMRI/'
data_dir = resdir + 'raw_data/'
output_dir = resdir + 'out/'
tmpdir = resdir + 'tmp/'
template = op.join(resdir,'template','MNI152_T1_2mm_brain.nii.gz')
template_3mm = op.join(resdir,'template','MNI152_T1_3mm_brain.nii.gz')
anatomy = 'FFE' # or TFE
functional_data = 'FS1' # other options: FS2, RS1, RS2 (four data-sets per subject; two fixed-state, two resting-state)
TR = 2.
fwhm = 5
highpass_sigma = 100/(2*TR)
lowpass_sigma = 12.5/(2*TR)

# next, define a subject list (each a directory)
subject_list = []
#subject_list = ['n3806','n3810','n3839','n3859','n3860','n3875','n3876','n3941','n3942','n3964',
#		'n3965','n3978','n3979','n3989','n3990','n4024','n4025','n4077','n4092','n4093',
#		'n4121','n4122','n4130','n4131','n4143','n4144','n4166','n4167','n4171','n4173',
#		'n4190','n4204','n4225','n4238','n4239','n4253','n4254','n4268','n4281','n4300',
#		'n4301','n4322','n4323','n4334','n4335','n4349','n4383','n4384','n4412','n4470',
#		'n4482','n4483','n4501','n4502']
# Map field names to individual subject runs: we have both functional data as anatomical data.

info = dict(func=[['subject_id', [functional_data]]],
            struct=[['subject_id',[anatomy]]])
#
"""
Set up parameters for the resting state preprocessing workflow.
"""
restingflow = create_resting_preproc()
restingflow.inputs.inputspec.num_noise_components = 6
restingflow.inputs.inputspec.highpass_sigma = highpass_sigma
restingflow.inputs.inputspec.lowpass_sigma = lowpass_sigma
restingflow.inputs.inputspec.template = template
restingflow.inputs.inputspec.template_3mm = template_3mm
infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),
                     name="infosource")
infosource.iterables = ('subject_id', subject_list)
datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=['func', 'struct']),
                     name = 'datasource')
datasource.inputs.base_directory = data_dir
datasource.inputs.template = '%s/%s.nii*'
datasource.inputs.template_args = info
datasource.inputs.sort_filelist = True

datasink = pe.Node(interface=nio.DataSink(parameterization=False),
                   name='datasink')
datasink.inputs.base_directory = output_dir

"""
Set up complete workflow
------------------------
"""
l1pipeline = pe.Workflow(name= "resting")
l1pipeline.base_dir = tmpdir
l1pipeline.connect([(infosource, datasource, [('subject_id', 'subject_id')]),
                    (datasource, restingflow, [('func', 'inputspec.func'),
                                               ('struct','inputspec.struct')]),
                    (infosource, datasink, [('subject_id', 'container'),
                                            (('subject_id', get_substitutions),
                                              'substitutions')]),
                    (restingflow, datasink, [('outputspec.noise_mask_file',
                                              '@noisefile'),
                                              ('outputspec.filtered_file',
                                               '@filteredfile')
                                              ])
               ])
l1pipeline.write_graph(graph2use='colored')
# l1pipeline.run(plugin='MultiProc', plugin_args={'n_procs' : 10})
# l1pipeline.run()
