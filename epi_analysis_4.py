# Step by step preparation for epilepsy fMRI analysis
from pathlib import Path
import os

# 0. SETUP
# List of patients to be analysed and the information belonging to them
patients = [
    {
        'patient_code': '',                     # Name code consisting of first two letters of first name and first two letters of last name
        'last_name': '',                     # Last nime as written in file path
        'first_name': '',                   # First name as written in file path
        'date': '',                           # Date as written in file path
        'paradigms': ['AMY_short', 'WFA'],         # List of paradigms that the patient has completed
        }, 
            ]

# Choose which functions should run
rename_dicoms   = False
create_vmr      = False
create_fmr      = True
coregister      = True
create_vtcs     = True

for patient in patients:
    # Set variables based on dict
    patient_code    = patient['patient_code']
    last_name       = patient['last_name']
    first_name      = patient['first_name']
    date            = patient['date']
    fmr_files       = patient['paradigms']

    patient_name    = f'{last_name}^{first_name}'
    datapath        = f'C:/BV_test/{last_name}_{first_name}_{date}/{first_name}_{last_name}_{date}_RENAME/'

    def get_file_path(fmr_type):
        """Checks if the file exists or if we need to add a space after the patient name"""
        new_path    = f'{datapath}{patient_code}_{fmr_type}'
        first_file = None
        for currentdir, subdirs, filenames in os.walk(new_path):
            for filename in filenames:
                if filename.endswith('001-0001-00001.dcm'):
                    first_file = os.path.join(currentdir, filename)
                    print(first_file, Path(first_file).exists())
                    return first_file
        if first_file == None:
            print(f'First file not found in path {new_path}')

    def select_fmr_vals(type):
        """
        Sets the correct values for each fmr and outputs all necessary values as a dictonary. 
        If new fMRI paradigms/protocols are introduced, their values must be added to this function. 
        """
        if type == 'AMY_short':
            protocol    = 'amygdala'
            volumes     = 272
            skip        = 16
            slices      = 24

        elif type == 'SPRACHE_SEM' or type == 'WFA' or type == 'SPRACHE_PHON':
            protocol = 'language'
            volumes     = 110
            skip        = 10
            slices      = 32
        
        elif type == 'HT':
            protocol = 'memory'
            volumes     = 210
            skip        = 10
            slices      = 32
        else:
            print(f'Error: Unknown fmr paradigm "{fmr}"')

        protocol_path = f'C:/BV_test/03_Stimulationsprotokolle/{protocol}.prt'

        fmr_vals = {'protocol': protocol_path, 'volumes': volumes, 'skip': skip, 'slices': slices}
        print(fmr_vals)
        return fmr_vals

    # 0.5 RENAME DICOM FILES
    # Dicom files are renamed according to BVs standards
    if rename_dicoms:
        folders = ['MPR'] + fmr_files
        for i in folders:
            print(f'Renaming dicoms for {i}')
            print(f'{datapath}{patient_code}_{i}')
            brainvoyager.rename_dicoms(f'{datapath}{patient_code}_{i}')

    # 1. Create VMR
    if create_vmr:
        first_file = get_file_path('MPR')
        print(f'Creating VMR using {first_file}')

        docVMR = brainvoyager.create_vmr('DICOM', first_file, 176, False, 256, 256, 2)

        # 1.5 VMR PREPROCESSING
        # Intensity inhomogeinities correction
        docVMR.correct_intensity_inhomogeneities_ext(False, 2, 0.25, 0.3, 3)
        print('IIHC')

        # Transform to isovoxel
        docVMR.transform_to_std_isovoxel(2, f'{patient_code}_MPR_ISO.vmr')
        docVMR_ISO = brainvoyager.open_document(f'{datapath}{patient_code}_MPR/{patient_code}_MPR_ISO.vmr')

        # Transform to sag -> only happens if needed (which it usually isn't with our files)
        sag_orientation = docVMR_ISO.transform_to_std_sag(f'{patient_code}_MPR_ISO_SAG.vmr')
        
    # 2. CREATE FMRS
    if create_fmr:
        # Loop through the completed paradigms
        for fmr in fmr_files:
            fmr_vals = select_fmr_vals(fmr)
            first_file = get_file_path(fmr)

            # fmr_dicom = brainvoyager.create_fmr_dicom(first_file, f'{patient_code}_{fmr}', f'{datapath}{patient_code}_{fmr}', fmr_vals[0])
            fmr_dicom = brainvoyager.create_fmr('DICOM', first_file, fmr_vals['volumes'], fmr_vals['skip'], True, fmr_vals['slices'], f'{patient_code}_{fmr}', False, 1,1,1, f'{datapath}{patient_code}_{fmr}')
            print(fmr_dicom)

            # 2.5 FMR Preprocessing
            # Spatial smoothing (AMY only)
            if fmr == 'AMY_short':
                fmr_dicom.smooth_spatial(4, 'mm')

            # Slice time correction
            fmr_dicom.correct_slicetiming(1, 2, 1)

            # 3D Motion correction
            fmr_dicom.correct_motion_ext(0, 2, False, 100, False, True)

            # Temporal Filtering 
            fmr_dicom.filter_temporal_highpass_fft(3, 1)

            # Link protocol
            link_protocol_ok = fmr_dicom.link_protocol(fmr_vals['protocol'])
            fmr_dicom.save()
            if link_protocol_ok:
                print('successfully linked protocol')
            else:
                print('failed to link protocol')

    # 3. Coregister VMR-FMRs
    if coregister:
        base_vmr = brainvoyager.open_document(f'{datapath}{patient_code}_MPR/{patient_code}_MPR_ISO_SAG.vmr')
        if base_vmr == None:
            base_vmr = brainvoyager.open_document(f'{datapath}{patient_code}_MPR/{patient_code}_MPR_ISO.vmr')

        for fmr in fmr_files:
            print(fmr)
            vmr = base_vmr
            fmr_doc = f'{datapath}{patient_code}_{fmr}/{patient_code}_{fmr}.fmr'
            coreg = base_vmr.coregister_fmr_to_vmr(fmr_doc, False, False)
            print(coreg)

    # 4. Create VTCs
    if create_vtcs:
        resolution      = 2
        interpolation   = 2  # 0 -> nearest neighbour, 1 -> trilinear, 2 -> sinc
        exttal          = False

        for fmr in fmr_files:
            fmr_path = f'{datapath}{patient_code}_{fmr}'
            for currentdir, subdirs, filenames in os.walk(fmr_path):
                vtcfiles = dict() # {'VMR', 'FMR', 'IA', 'FA'} 
                brainvoyager.print_to_log('Reset files...')
                brainvoyager.print_to_log('Current directory: ' + currentdir)
                for filename in filenames:           
                    print('Current directory: ' + currentdir) 
                    if filename.endswith('.fmr') and not filename.endswith('_firstvol.fmr'): 
                        vtcfiles['FMR'] = os.path.join(currentdir, filename)
                        brainvoyager.print_to_log('FMR found: ' + vtcfiles['FMR'])
                    if filename.endswith('.vmr'):     
                        vtcfiles['VMR'] = os.path.join(currentdir, filename)
                        brainvoyager.print_to_log('VMR found: ' + vtcfiles['VMR'])
                    if filename.endswith('_IA.trf'):     
                        vtcfiles['IA'] = os.path.join(currentdir, filename)     
                        brainvoyager.print_to_log('Initial alignment file found: ' + vtcfiles['IA'])       
                    if filename.endswith('_FA.trf'):     
                        vtcfiles['FA'] = os.path.join(currentdir, filename)
                        brainvoyager.print_to_log('Fine alignment file found: ' + vtcfiles['FA']) 
                    if (len(vtcfiles.values()) == 4):    
                        brainvoyager.print_to_log('Number of unique file types found: ' + str(len(vtcfiles)))    
                        brainvoyager.print_to_log('All files present, proceed with creating VTC...')
                        doc = brainvoyager.open_document(vtcfiles['VMR'])
                        newname = os.path.join(currentdir, vtcfiles['FMR'].split(".fmr", 1)[0]+'_res'+str(resolution)+'_interp'+str(interpolation)+"_native.vtc")
                        doc.vtc_creation_extended_tal_space = exttal
                        brainvoyager.print_to_log('Creating ' + newname + ' from:' +'\n' +'\t'+vtcfiles['FMR'] + '\n' +'\t'+ vtcfiles['IA'] + '\n'+ '\t'+ vtcfiles['FA'])
                        ok = doc.create_vtc_in_native_space(vtcfiles['FMR'], vtcfiles['IA'], vtcfiles['FA'], newname, resolution, interpolation)
                        # doc.close() 
                        vtcfiles = dict() # reset   
