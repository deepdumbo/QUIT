#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities.

Currently implemented:
    - qidespot1

Requires that the QUIT tools are in your your system path

"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
import json

############################ qidespot1 ############################

class QIDespot1InputSpec(CommandLineInputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to SPGR data')

    param_file = File(desc='Parameter .json file', position=1, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=1, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    
class QIDespot1OutputSpec(TraitedSpec):
    t1_map = File(desc="Path to T1 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QIDespot1(CommandLine):
    """
    Run DESPOT1 analysis with qidespot1

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1
    >>> d1 = QIDespot1(prefix='nipype_', param_file='spgr_params.json')
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]} }
    >>> d1 = QIDespot1(prefix='nipype_', param_dict=params)
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    """

    _cmd = 'qidespot1'
    input_spec = QIDespot1InputSpec
    output_spec = QIDespot1OutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot1, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t1_map'] = prefix + 'D1_T1.nii.gz'
        outputs['pd_map'] = prefix + 'D1_PD.nii.gz'
        
        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'D1_residual.nii.gz'
        
        return outputs

############################ qidespot1sim ############################

class QIDespot1SimInputSpec(CommandLineInputSpec):
    # Inputs
    spgr_file = File(exists=False, argstr='%s', mandatory=True,
        position=0, desc='Path to write SPGR/FLASH image')

    param_file = File(desc='Parameter .json file', position=1, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=1, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    noise = traits.Float(desc='Noise level to add to simulation', argstr='--simulate=%f', default_value=0.0, usedefault=True)
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    
class QIDespot1SimOutputSpec(TraitedSpec):
    spgr_image = File(desc="Path to SPGR/FLASH image")

class QIDespot1Sim(CommandLine):
    """
    Run DESPOT1 simulation with qidespot1

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]},
                  'T1File': 'D1_T1.nii.gz',
                  'PDFile': 'D1_PD.nii.gz'}
    >>> d1sim = QIDespot1Sim(prefix='nipype_', param_dict=params)
    >>> d1sim.inputs.spgr_file = 'SPGR.nii.gz'
    >>> d1sim_res = d1.run()
    >>> print(d1sim_res.outputs)

    """

    _cmd = 'qidespot1'
    input_spec = QIDespot1SimInputSpec
    output_spec = QIDespot1SimOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot1Sim, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['spgr_image'] = self.inputs.spgr_file
        return outputs

############################ qidespot1hifi ############################

class QIDespot1HifiInputSpec(CommandLineInputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to SPGR data')

    mprage_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to MPRAGE data')

    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QIDespot1HifiOutputSpec(TraitedSpec):
    t1_map = File(desc="Path to T1 map")
    pd_map = File(desc="Path to PD map")
    b1_map = File(desc="Path to B1 map")
    residual_map = File(desc="Path to residual map")

class QIDespot1Hifi(CommandLine):
    """
    Calculate T1 & B1 map with the DESPOT1-HIFI method

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1Hifi
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]}, 
                  'MPRAGE': { 'FA': 5, 'TR': 5E-3, 'TI': 0.45, 'TD': 0, 'eta': 1, 'ETL': 64, 'k0': 0 }}
    >>> hifi = QIDespot1Hifi(prefix='nipype_', param_dict=params)
    >>> hifi.inputs.spgr_file = 'SPGR.nii.gz'
    >>> hifi.inputs.mprage_file = 'MPRAGE.nii.gz'
    >>> hifi_res = hifi.run()
    >>> print(hifi_res.outputs)

    """

    _cmd = 'qidespot1hifi'
    input_spec = QIDespot1HifiInputSpec
    output_spec = QIDespot1HifiOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot1Hifi, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t1_map'] = prefix + 'HIFI_T1.nii.gz'
        outputs['pd_map'] = prefix + 'HIFI_PD.nii.gz'
        outputs['b1_map'] = prefix + 'HIFI_B1.nii.gz'

        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'HIFI_residual.nii.gz'
        
        return outputs

############################ qidespot1hifisim ############################

class QIDespot1HifiSimInputSpec(CommandLineInputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path for output SPGR image')

    mprage_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path for output MPRAGE image')

    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    noise = traits.Float(desc='Noise level to add to simulation', argstr='--simulate=%f', default_value=0.0, usedefault=True)
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')

    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QIDespot1HifiSimOutputSpec(TraitedSpec):
    spgr_image = File(desc="Output SPGR image")
    mprage_image = File(desc="Output MPRAGE image")

class QIDespot1HifiSim(CommandLine):
    """
    Simulate SPGR/FLASH and MPRAGE images using DESPOT1-HIFI model

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1HifiSim
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]}, 
                  'MPRAGE': { 'FA': 5, 'TR': 5E-3, 'TI': 0.45, 'TD': 0, 'eta': 1, 'ETL': 64, 'k0': 0 },
                  'PDFile': 'PD.nii.gz',
                  'T1File': 'T1.nii.gz'}
    >>> hifisim = QIDespot1HifiSim(prefix='nipype_', param_dict=params)
    >>> hifisim.inputs.spgr_file = 'SPGR.nii.gz'
    >>> hifisim.inputs.mprage_file = 'MPRAGE.nii.gz'
    >>> hifisim_res = hifi.run()
    >>> print(hifisim_res.outputs)

    """

    _cmd = 'qidespot1hifi'
    input_spec = QIDespot1HifiSimInputSpec
    output_spec = QIDespot1HifiSimOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot1HifiSim, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['spgr_image'] = self.inputs.spgr_file
        outputs['mprage_image'] = self.inputs.mprage_file
        return outputs

############################ qidespot2 ############################

class QIDespot2InputSpec(CommandLineInputSpec):
    # Inputs

    t1_map = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to T1 map')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
        position=1, desc='Path to SSFP data')
    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)
    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    ellipse = traits.Bool(desc="Data is ellipse geometric solution", argstr='-e')
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T2 = traits.Float(desc='Clamp T2 between 0 and value', argstr='--clampT1=%f')
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QIDespot2OutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QIDespot2(CommandLine):
    """
    Run DESPOT2 analysis with qidespot2

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2
    >>> d1 = QIDespot2(prefix='nipype_', param_file='ssfp_params.json')
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]} }
    >>> d2 = QIDespot2(prefix='nipype_', param_dict=params)
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    """

    _cmd = 'qidespot2'
    input_spec = QIDespot2InputSpec
    output_spec = QIDespot2OutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot2, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t2_map'] = prefix + 'D2_T2.nii.gz'
        outputs['pd_map'] = prefix + 'D2_PD.nii.gz'
        
        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'D2_residual.nii.gz'
        
        return outputs

############################ qidespot2sim ############################

class QIDespot2SimInputSpec(CommandLineInputSpec):
    # Inputs

    t1_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path for input T1 map')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
        position=1, desc='Path for output SSFP data')
    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)
    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    noise = traits.Float(desc='Noise level to add to simulation', argstr='--simulate=%f', default_value=0.0, usedefault=True)
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    ellipse = traits.Bool(desc="Data is ellipse geometric solution", argstr='-e')
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QIDespot2SimOutputSpec(TraitedSpec):
    ssfp_image = File(desc="Path to SSFP image")

class QIDespot2Sim(CommandLine):
    """
    Run DESPOT2 simulation

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2Sim
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]},
                  'PDFile': 'PD.nii.gz',
                  'T1File': 'T1.nii.gz' }
    >>> d2sim = QIDespot2Sim(prefix='nipype_', param_dict=params)
    >>> d2sim.inputs.in_file = 'SSFP.nii.gz'
    >>> d2sim.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2sim_res = d2.run()
    >>> print(d2sim_res.outputs)

    """

    _cmd = 'qidespot2'
    input_spec = QIDespot2SimInputSpec
    output_spec = QIDespot2SimOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot2Sim, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['ssfp_image'] = self.inputs.ssfp_file
        return outputs

############################ qidespot2fm ############################

class QIDespot2FMInputSpec(CommandLineInputSpec):
    # Inputs

    t1_map = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to T1 map')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
        position=1, desc='Path to SSFP data')
    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)
    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    asym = traits.Bool(desc="Fit asymmetric (+/-) off-resonance frequency", argstr='-A')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QIDespot2FMOutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    f0_map = File(desc="Path to f0 (off-resonance) map")
    residual_map = File(desc="Path to residual map")

class QIDespot2FM(CommandLine):
    """
    Run DESPOT2-FM analysis

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2FM
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,10,50,50], 'PhaseInc':[180,180,0,0] }
    >>> fm = QIDespot2FM(prefix='nipype_', param_dict=params)
    >>> fm.inputs.in_file = 'SSFP.nii.gz'
    >>> fm.inputs.t1_file = 'D1_T1.nii.gz'
    >>> fm_res = fm.run()
    >>> print(fm_res.outputs)

    """

    _cmd = 'qidespot2fm'
    input_spec = QIDespot2FMInputSpec
    output_spec = QIDespot2FMOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot2FM, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t2_map'] = prefix + 'FM_T2.nii.gz'
        outputs['pd_map'] = prefix + 'FM_PD.nii.gz'
        outputs['f0_map'] = prefix + 'FM_f0.nii.gz'
        
        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'FM_residual.nii.gz'
        
        return outputs

############################ qimcdespot ############################

class QIDespot2FMInputSpec(CommandLineInputSpec):
    # Inputs

    t1_map = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to T1 map')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
        position=1, desc='Path to SSFP data')
    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)
    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    asym = traits.Bool(desc="Fit asymmetric (+/-) off-resonance frequency", argstr='-A')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QIDespot2FMOutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    f0_map = File(desc="Path to f0 (off-resonance) map")
    residual_map = File(desc="Path to residual map")

class QIDespot2FM(CommandLine):
    """
    Run DESPOT2-FM analysis

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2FM
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,10,50,50], 'PhaseInc':[180,180,0,0] }
    >>> fm = QIDespot2FM(prefix='nipype_', param_dict=params)
    >>> fm.inputs.in_file = 'SSFP.nii.gz'
    >>> fm.inputs.t1_file = 'D1_T1.nii.gz'
    >>> fm_res = fm.run()
    >>> print(fm_res.outputs)

    """

    _cmd = 'qidespot2fm'
    input_spec = QIDespot2FMInputSpec
    output_spec = QIDespot2FMOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QIDespot2FM, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t2_map'] = prefix + 'FM_T2.nii.gz'
        outputs['pd_map'] = prefix + 'FM_PD.nii.gz'
        outputs['f0_map'] = prefix + 'FM_f0.nii.gz'
        
        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'FM_residual.nii.gz'
        
        return outputs