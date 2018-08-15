%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.cfg_mkdir.parent = {[pwd '/']};
matlabbatch{1}.cfg_basicio.cfg_mkdir.name = 'seg_result';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1) = cfg_dep;
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).tname = 'Directory';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).tgt_spec{1}(1).value = 'dir';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).sname = 'Make Directory: Make Directory ''seg_result''';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).src_output = substruct('.','dir');
matlabbatch{3}.spm.util.dicom.data = '<UNDEFINED>';
matlabbatch{3}.spm.util.dicom.root = 'flat';
matlabbatch{3}.spm.util.dicom.outdir(1) = cfg_dep;
matlabbatch{3}.spm.util.dicom.outdir(1).tname = 'Output directory';
matlabbatch{3}.spm.util.dicom.outdir(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{3}.spm.util.dicom.outdir(1).tgt_spec{1}(1).value = 'dir';
matlabbatch{3}.spm.util.dicom.outdir(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.util.dicom.outdir(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.util.dicom.outdir(1).sname = 'Make Directory: Make Directory ''seg_result''';
matlabbatch{3}.spm.util.dicom.outdir(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.util.dicom.outdir(1).src_output = substruct('.','dir');
matlabbatch{3}.spm.util.dicom.convopts.format = 'nii';
matlabbatch{3}.spm.util.dicom.convopts.icedims = 0;
matlabbatch{4}.spm.spatial.preproc.data(1) = cfg_dep;
matlabbatch{4}.spm.spatial.preproc.data(1).tname = 'Data';
matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(1).value = 'image';
matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.preproc.data(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.preproc.data(1).sname = 'DICOM Import: Converted Images';
matlabbatch{4}.spm.spatial.preproc.data(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.preproc.data(1).src_output = substruct('.','files');
matlabbatch{4}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{4}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{4}.spm.spatial.preproc.output.CSF = [0 0 1];
matlabbatch{4}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{4}.spm.spatial.preproc.output.cleanup = 2;
matlabbatch{4}.spm.spatial.preproc.opts.tpm = {
                                               'C:\Users\LAB-690\Desktop\104project\spm8\tpm\grey.nii'
                                               'C:\Users\LAB-690\Desktop\104project\spm8\tpm\white.nii'
                                               'C:\Users\LAB-690\Desktop\104project\spm8\tpm\csf.nii'
                                               };
matlabbatch{4}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{4}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{4}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{4}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{4}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{4}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{4}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{4}.spm.spatial.preproc.opts.msk = {''};
