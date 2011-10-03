"""
Author: Robert Louis Cloud, rcloud@gmail.com
Modified: 11/08/2011

script.py
"""

import os
import nipype.pipeline.engine as pe
from nipype.interfaces import fsl
import nibabel
import sys
import datetime

def py_erode(in_file, delta_band):
    import nipype.interfaces.matlab as matlab
    import os
    mlab = matlab.MatlabCommand()
#    erode_val = 10
    mlab.inputs.script = "erode_mask_cereb_flr('%s', '%s')" % (in_file, delta_band)
    out = mlab.run()
    return os.path.join(os.getcwd(), 'out_mask.nii')

def usage():
    print """
             Usage: python script.py subject_number upper_threshold lower_threshold delta_band 
	     
	     example: python script.py c011v0 2.0 2.5 15

	     2.0 # standard dev. from peripheral mask
	     2.5 # lower threshold
	     15 is the number of voxels removed from the outer edge of cortical mask	
	"""
    sys.exit(-1)
            
if len(sys.argv) != 5:
    usage()
subject_num = sys.argv[1]
zval30 = float(sys.argv[2])
zval35 = float(sys.argv[3])
delta_band = int(sys.argv[4])


cur_dir = os.getcwd()
cerebellar_mask = cur_dir + "/masks/cerebellar_mask.nii"
brainstem_mask  = cur_dir + "/masks/brainstem_mask.nii"
periph_mask     = "/usr/local/fsl/data/standard/MNI152_T1_2mm_strucseg_periph.nii.gz"
mprage_path     = cur_dir + "/flair/" + subject_num + "/mprage.nii"
flair_path      = cur_dir + "/flair/" + subject_num + "/flair.nii"
mrn_path        = cur_dir + "/flair/" + subject_num + "/mrn/"

fsl.FSLCommand.set_default_output_type('NIFTI')

#first setup logging
filename = "/%s_log.txt" %subject_num
log = open(cur_dir + filename, "a")
log.write("****\n")
now = datetime.datetime.now()
log.write(now.strftime("%Y-%m-%d %H:%M"))
log.write("\nlower threshold:\t " + sys.argv[2])
log.write("\nupper threshold:\t " + sys.argv[3])
log.write("\ndelta band:\t\t " + sys.argv[4])
log.write("\n")




"""
1) create a workflow to contain the nodes
files will be created in this subdirectory of the current directory
"""
wf_name = "%s_workflow" % subject_num
of_name = "%s_output" % subject_num
workflow = pe.Workflow(name= wf_name)
workflow.base_dir = '.'

"""
2)Create the nodes to make the cerebellar and brainstem masks for the subject
"""
subject_cerebellar_mask = pe.Node(interface = fsl.FLIRT(in_file=cerebellar_mask,
                                                reference=mprage_path,
                                                in_matrix_file = mrn_path + "I2std_inv.mat",
                                                apply_xfm = True),
                          name = 'subject_cerebellar_mask')

subject_brainstem_mask = pe.Node(interface = fsl.FLIRT(in_file=brainstem_mask,
                                                reference=mprage_path,
                                                in_matrix_file = mrn_path + "I2std_inv.mat",
                                                apply_xfm = True),
                          name = 'subject_brainstem_mask')


""" 
3)Create the nodes for inverted subject mask for the cerebellum and brainstem 
  and connect them to the output of the flirt created cerebellum and brainstem masks.  
creating inverse transformation matrices
"""

subject_cerebellar_inv_mask = pe.Node(interface=fsl.ImageMaths(op_string='-bin -mul -1 -add 1'),                                                  
                          name = 'subject_cerebellar_inv_mask')
subject_brainstem_inv_mask = pe.Node(interface=fsl.ImageMaths(op_string='-bin -mul -1 -add 1'),                                                 
                         name = 'subject_brainstem_inv_mask')

workflow.connect(subject_cerebellar_mask, 'out_file', subject_cerebellar_inv_mask, 'in_file')
workflow.connect(subject_brainstem_mask, 'out_file', subject_brainstem_inv_mask, 'in_file')

"""
4) perform thresholding to create subject_vent_inv_mask
where did this -thr .3 come from?
"""
subject_vent_inv_mask = pe.Node(interface=fsl.ImageMaths(in_file= mrn_path + "I_stdmaskbrain_pve_0_segvent.nii",
                                                         op_string="-thr .3 -add 1 -uthr 1 -bin"),
                            name = 'subject_vent_inv_mask')





"""
5) Apply masks to create subject cerebral mask
input: I_stdmaskbrain
mask1: subject_cerebellar_inv_mask
mask2: subject_brainstem_inv_mask
mask3: subject_vent_inv_mask
option: -bin
output: subject_cerebral_mask
"""

apply_cerebellar_mask = pe.Node(interface = fsl.maths.ApplyMask(in_file = mrn_path + "I_stdmaskbrain.nii"),
                                name='apply_cerebellar_mask')

#connect the cerebellar_inv_mask to the mask_file input
workflow.connect(subject_cerebellar_inv_mask, 'out_file', apply_cerebellar_mask, 'mask_file')


apply_brainstem_mask = pe.Node(interface = fsl.maths.ApplyMask(),
                               name = 'apply_brainstem_mask')
#connect the output of apply_cerebellar_mask to the input of apply_brainstem_mask
workflow.connect(apply_cerebellar_mask, 'out_file', apply_brainstem_mask, 'in_file')

#connect the inverse brainstem mask to the mask file of apply_brainstem mask
workflow.connect(subject_brainstem_inv_mask, 'out_file', apply_brainstem_mask, 'mask_file')


apply_vent_mask = pe.Node(interface = fsl.maths.ApplyMask(),
                          name = 'apply_vent_mask')
#connect the output of apply brainstem to the input of apply vent
workflow.connect(apply_brainstem_mask, 'out_file', apply_vent_mask, 'in_file')

#connect the output of vent_inv_mask to the mask file of apply vent
workflow.connect(subject_vent_inv_mask, 'out_file', apply_vent_mask, 'mask_file')

#binarize to get final cerebral mask output
subject_cerebral_mask = pe.Node(interface = fsl.maths.UnaryMaths(operation='bin'),
                                name='subject_cerebral_mask')
workflow.connect(apply_vent_mask, 'out_file', subject_cerebral_mask, 'in_file')


"""
6) Perform another masking with the just created subject_cerebral_mask
   to produce subject_cerebral_pve_2
"""
subject_cerebral_pve_2 = pe.Node(interface = fsl.maths.ApplyMask(in_file = mrn_path + "I_stdmaskbrain_pve_2.nii"),                                                           
                                 name = 'subject_cerebral_pve_2')

workflow.connect(subject_cerebral_mask, 'out_file', subject_cerebral_pve_2, 'mask_file')                                                          


                                                             
"""
7) Thresholding to create cerebral tissue mask
operations: -thr .1 -add 1 -uthr
mask with subject_cerebral_mask
The probability map for CSF output by SIENAX was used to identify the voxels having at least a 10% chance of being CSF in the MPRAGE. 
These voxels were removed from cerebral brain mask generating a cerebral tissue mask for the MPRAGE image.
"""
cerebral_tissue_ops = pe.Node(interface=fsl.ImageMaths(in_file = mrn_path + "I_stdmaskbrain_pve_0.nii",
                                                       op_string = "-thr .1 -add 1 -uthr 1 -bin "),
                              name='cerebral_tissue_ops')
workflow.add_nodes([cerebral_tissue_ops])

subject_cerebral_tissue_mask = pe.Node(interface=fsl.maths.ApplyMask(), 
                               name='subject_cerebral_tissue_mask')
workflow.connect(cerebral_tissue_ops, 'out_file', subject_cerebral_tissue_mask, 'in_file')
workflow.connect(subject_cerebral_mask, 'out_file', subject_cerebral_tissue_mask, 'mask_file')


"""
8) Generate subject_cerebral_WM_mask
"""
subject_cerebral_WM_mask = pe.Node(interface=fsl.ImageMaths(op_string = "-thr .5 -bin"),
                                   name='subject_cerebral_WM_mask')
workflow.connect(subject_cerebral_pve_2, 'out_file', subject_cerebral_WM_mask, 'in_file')


"""
9) Generate a mask for peripheral cortex in MPRAGE space
"""
subject_periph_mask = pe.Node(interface = fsl.FLIRT(in_file=periph_mask,
                                                reference=mprage_path,
                                                in_matrix_file = mrn_path + "I2std_inv.mat",
                                                apply_xfm = True),
                          name = 'subject_periph_mask')

subject_periph_mask_binarized = pe.Node(interface=fsl.maths.UnaryMaths(operation='bin'),
                          name = 'subject_periph_mask_binarized')
workflow.connect(subject_periph_mask, 'out_file', subject_periph_mask_binarized, 'in_file')


"""
10) Generate WM PVE image of peripheral cortex region
"""

#mask subject_cerebral_pve_2 with subject_periph_mask_binarized
mask_cerebral_pve_2_with_periph = pe.Node(interface = fsl.maths.ApplyMask(),
                             name = 'mask_cerebral_pve_2_with_periph')
workflow.connect(subject_cerebral_pve_2, 'out_file', mask_cerebral_pve_2_with_periph, 'in_file')
workflow.connect(subject_periph_mask_binarized, 'out_file', mask_cerebral_pve_2_with_periph, 'mask_file')

#create the periph WM mask by binarizing
subject_cerebral_periph_WM_mask = pe.Node(interface = fsl.maths.UnaryMaths(operation='bin'),
                                          name = 'subject_cerebral_periph_WM_mask')
workflow.connect(mask_cerebral_pve_2_with_periph, 'out_file', subject_cerebral_periph_WM_mask, 'in_file')

                 

"""
11) invert subject_cerebral_WM_mask to create subject_cerebral_WM_inv_mask
"""
subject_cerebral_WM_inv_mask = pe.Node(interface=fsl.ImageMaths(op_string = "-sub 1 -mul -1 -bin"),
                                       name='subject_cerebral_WM_inv_mask')
workflow.connect(subject_cerebral_WM_mask, 'out_file', subject_cerebral_WM_inv_mask, 'in_file')

"""
12) mask subject_cerebral_periph_WM_mask with subject_cerebral_WM_inv_mask
    to produce subject_cerebral_periph_mixed_mask
"""


subject_cerebral_periph_mixed_mask = pe.Node(interface=fsl.maths.ApplyMask(),
                                             name='subject_cerebral_periph_mixed_mask')
workflow.connect(subject_cerebral_periph_WM_mask, 'out_file', subject_cerebral_periph_mixed_mask, 'in_file')
workflow.connect(subject_cerebral_WM_inv_mask, 'out_file', subject_cerebral_periph_mixed_mask, 'mask_file')



"""
#13) Register the mprage cerebral mask, WM mask, and tissue mask to the FLAIR image
# this is the most time consuming step after sienax
This generates a nifti image of the mprage image to flair space as well as a transformation 
matrix which can be used to transform the MPRAGE generated masks to FLAIR space 
and then applied to the FLAIR image.
"""

register_mprage = pe.Node(interface=fsl.FLIRT(in_file=mprage_path,
                                              reference=flair_path,                                         
                                              bins=256,
                                              cost='normmi',
                                              searchr_x=[-5, 5],
                                              searchr_y=[-5, 5],
                                              searchr_z=[-5, 5],
                                              dof=12,
                                              interp='sinc',
                                              sinc_width=7,
                                              sinc_window='hanning'),
                          name='register_mprage')
                                              
workflow.add_nodes([register_mprage])                                              


"""
For the following, the steps don't take long and we can manipulate the 
threshold values
"""

#14) Register the cerebral mask


#create flair cerebral mask
flair_cerebral_mask = pe.Node(interface = fsl.FLIRT(reference=flair_path,
                                                apply_xfm = True),
                          name = 'flair_cerebral_mask')
workflow.connect(subject_cerebral_mask, 'out_file', flair_cerebral_mask, 'in_file')
workflow.connect(register_mprage, 'out_matrix_file', flair_cerebral_mask, 'in_matrix_file')

#we have two masks because the second has the threshold of 50% applied
#this is a binary threshold where only the values greater than .5 are included.
flair_cerebral_mask_deux = pe.Node(interface=fsl.ImageMaths(op_string = "-thr .5 -bin"),
                                   name='flair_cerebral_mask_deux')
#workflow.connect(flair_cerebral_mask, 'out_file', flair_cerebral_mask_deux, 'in_file')
workflow.connect([(flair_cerebral_mask, flair_cerebral_mask_deux,[(('out_file'), 'in_file')]),])
#take flair_cerebral_mask_deux and mask out the edges to remove WMH artifacts


"""
matlab_erode = pe.Node(interface=erode.ConmapTxt2Mat(script=matlab_script),
                name = 'matlab_erode')

workflow.connect(flair_cerebral_mask_deux, 'out_file', matlab_erode, 'in_file')
"""


#apply flair cerebral mask to flair
flair_cerebral = pe.Node(interface = fsl.maths.ApplyMask(in_file = flair_path),
                     name = 'flair_cerebral')
workflow.connect([(flair_cerebral_mask_deux, flair_cerebral,[(('out_file', py_erode, delta_band),'mask_file')]),])
#workflow.connect(flair_cerebral_mask_deux, 'out_file', flair_cerebral, 'mask_file')


#15) Register the WM mask with flair

flair_cerebral_WM_mask = pe.Node(interface = fsl.FLIRT(reference=flair_path,
                                                apply_xfm = True),
                          name = 'flair_cerebral_WM_mask')

workflow.connect(subject_cerebral_WM_mask, 'out_file', flair_cerebral_WM_mask, 'in_file')
workflow.connect(register_mprage, 'out_matrix_file', flair_cerebral_WM_mask, 'in_matrix_file')

flair_cerebral_WM_mask_deux = pe.Node(interface=fsl.ImageMaths(op_string = "-thr .5 -bin"),
                                      name='flair_cerebral_WM_mask_deux')
workflow.connect(flair_cerebral_WM_mask, 'out_file', flair_cerebral_WM_mask_deux, 'in_file')

flair_cerebral_WM = pe.Node(interface = fsl.maths.ApplyMask(in_file = flair_path),
                            name='flair_cerebral_WM')
#workflow.connect(flair_cerebral_WM_mask_deux, 'out_file', flair_cerebral_WM, 'mask_file')
workflow.connect([(flair_cerebral_WM_mask_deux, flair_cerebral_WM,[(('out_file', py_erode, delta_band),'mask_file')]),])


#16) Register the cerebral tissue mask with flair

flair_cerebral_tissue_mask = pe.Node(interface = fsl.FLIRT(reference=flair_path,
                                                apply_xfm = True),
                          name = 'flair_cerebral_tissue_mask')
workflow.connect(subject_cerebral_tissue_mask, 'out_file', flair_cerebral_tissue_mask, 'in_file')
workflow.connect(register_mprage, 'out_matrix_file', flair_cerebral_tissue_mask, 'in_matrix_file')

flair_cerebral_tissue_mask_deux = pe.Node(interface=fsl.ImageMaths(op_string = "-thr .5 -bin"),
                                          name='flair_cerebral_tissue_mask_deux')
workflow.connect(flair_cerebral_tissue_mask, 'out_file', flair_cerebral_tissue_mask_deux, 'in_file')

flair_cerebral_tissue = pe.Node(interface = fsl.maths.ApplyMask(in_file = flair_path),
                                name='flair_cerebral_tissue')
workflow.connect([(flair_cerebral_tissue_mask_deux, flair_cerebral_tissue,[(('out_file', py_erode, delta_band),'mask_file')]),])
#workflow.connect(flair_cerebral_tissue_mask_deux, 'out_file', flair_cerebral_tissue, 'mask_file')



#17) Register the cerebral peripheral mixed mask
# likely we will want to change the thresholding values for this...


flair_cerebral_periph_mixed_mask = pe.Node(interface = fsl.FLIRT(reference=flair_path,
                                                apply_xfm = True),
                          name = 'flair_cerebral_periph_mixed_mask')
workflow.connect(subject_cerebral_periph_mixed_mask, 'out_file', flair_cerebral_periph_mixed_mask, 'in_file')
workflow.connect(register_mprage, 'out_matrix_file', flair_cerebral_periph_mixed_mask, 'in_matrix_file')

flair_cerebral_periph_mixed_mask_deux = pe.Node(interface=fsl.ImageMaths(op_string= "-thr .5 -bin"),
                                                name="flair_cerebral_periph_mixed_mask_deux")
workflow.connect(flair_cerebral_periph_mixed_mask, 'out_file', flair_cerebral_periph_mixed_mask_deux, 'in_file')

flair_cerebral_periph_mixed = pe.Node(interface = fsl.maths.ApplyMask(in_file = flair_path),
                                      name='flair_cerebral_periph_mixed')
workflow.connect(flair_cerebral_periph_mixed_mask_deux, 'out_file', flair_cerebral_periph_mixed, 'mask_file')











"""
18) Estimate WMH volume
  Calculate WMH threshold using peripheral cerebral mixed
  use nipype.interfaces.fsl.ImageStats
  we are using the flair_cerebral_periph_mixed to calculate the standard deviation
  
"""

mean_peri = pe.Node(interface=fsl.ImageStats(op_string = '-M'),
                    name = 'mean_peri')
workflow.connect(flair_cerebral_periph_mixed, 'out_file', mean_peri, 'in_file')

std_peri = pe.Node(interface=fsl.ImageStats(op_string = '-S'),
                   name = 'std_peri')
workflow.connect(flair_cerebral_periph_mixed, 'out_file', std_peri, 'in_file')


# we get back a graph containing all of the nodes with the output of their work
graph = workflow.run()
nodes = graph.nodes()
workflow.write_graph()

"""
this is really a nasty hack to do string comparison to get the index of the node in the
networkx data system, but it works
#get index of node mean_peri and std_peri
"""
mean_peri_string = str(mean_peri)
std_dev_string = str(std_peri)
flair_cerebral_tissue_string = str(flair_cerebral_tissue)
flair_cerebral_WM_string = str(flair_cerebral_WM)

"""
could just use a switch statement here 
"""

for i in range(len(nodes)):
    node_str = str(nodes[i])
    if node_str == mean_peri_string:
        mean_index = i
    elif node_str == std_dev_string:
        std_index = i
    elif node_str == flair_cerebral_tissue_string:
        tissue_index = i
    elif node_str == flair_cerebral_WM_string:
        WM_index = i
        
mean_peri_num = graph.nodes()[mean_index].result.outputs.out_stat
std_output_num = graph.nodes()[std_index].result.outputs.out_stat
WM_location = graph.nodes()[WM_index].result.outputs.out_file
tissue_location = graph.nodes()[tissue_index].result.outputs.out_file

log.write("\n" + str(mean_peri_num))
log.write("\n" + str(std_output_num))

#We will be wanting to manipulate the std values of the following two steps
#get the three standard deviations from the peripheral mask level
#level30z_peri = mean_peri_num + 1.5 * std_output_num
level30z_peri = mean_peri_num + zval30 * std_output_num

#get the 3.5 standard deviation from the peripheral mask level
#level35z_peri = mean_peri_num + 2.2 * std_output_num
level35z_peri = mean_peri_num + zval35 * std_output_num


#Apply threshold to data
#create a new workflow for the final processing steps
output = pe.Workflow(name= of_name)
output.base_dir = '.'


flair_tissue_peri_thr35z = pe.Node(interface=fsl.ImageMaths(in_file = tissue_location, op_string="-thr " + str(level35z_peri) + " -bin"),
                                   name='flair_tissue_peri_thr35z')

flair_WM_peri_thr30z = pe.Node(interface=fsl.ImageMaths(in_file = WM_location, op_string="-thr " + str(level30z_peri) + " -bin"),
                               name='flair_WM_peri_thr30z')

output.add_nodes([flair_tissue_peri_thr35z, flair_WM_peri_thr30z])

flair_WM_30plus35 = pe.Node(interface=fsl.maths.BinaryMaths(operation='add'), 
                            name='flair_WM_30plus35')

output.connect(flair_WM_peri_thr30z, 'out_file', flair_WM_30plus35, 'in_file')
output.connect(flair_tissue_peri_thr35z, 'out_file', flair_WM_30plus35, 'operand_file')

flair_WM_30plus35_binarized = pe.Node(interface=fsl.maths.UnaryMaths(operation='bin'),
                                      name='flair_WM_30plus35_binarized')
output.connect(flair_WM_30plus35, 'out_file', flair_WM_30plus35_binarized, 'in_file')

WMH_out_file = cur_dir + "/" + of_name + "/" + str(zval30) + "_" + str(zval35) + ".nii"
flair_WM_30plus35_binarized.inputs.out_file = WMH_out_file

output.run()
output.write_graph()

WMH_file = nibabel.load(WMH_out_file)
WMH_data = WMH_file.get_data()
WMH_voxels = len(WMH_data.nonzero()[0])
WMH_hdr = WMH_file.get_header()
dimensions = WMH_hdr.get_zooms()
WMH_volume = dimensions[0] * dimensions[1] * dimensions[2] * WMH_voxels
print "WMH Volume = " + str(WMH_volume) + " mm^3"
log.write("\nWMH Volume: \t\t" + str(WMH_volume) + " mm^3")
log.write("\nWMH Volume: \t\t" + str(WMH_volume / 1000) + " ml")
log.write("\n****")
log.write("\n\n")
log.close()
imvals = open(cur_dir + "/" + of_name + "/" + "volume_levels.txt", "a")
imvals.write("\n" + str(zval30) + "\t" + str(zval35) + "\t" + str(WMH_volume / 1000))
imvals.close()
