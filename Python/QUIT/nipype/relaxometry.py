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

class QiDespot1InputSpec(CommandLineInputSpec):
    
    # Inputs
    param_file = File(desc='Parameter .json file', position=-1, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=-1, 
        argstr='', mandatory=True, xor=['param_file'])
        
    in_file = File(exists=True, argstr='%s', mandatory=True,
        position=-2, desc='Path to SPGR data')
    
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
    
class QiDespot1OutputSpec(TraitedSpec):
    t1_map = File(desc="Path to T1 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QiDespot1(CommandLine):
    """
    Run DESPOT1 analysis with qidespot1

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot1
    >>> d1 = QiDespot1(prefix='nipype_', param_file='spgr_params.json')
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot1
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]} }
    >>> d1 = QiDespot1(prefix='nipype_', param_dict=params)
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    """

    _cmd = 'qidespot1'
    input_spec = QiDespot1InputSpec
    output_spec = QiDespot1OutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiDespot1, self)._format_arg(name, spec, value)

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

############################ qidespot2 ############################

class QiDespot2InputSpec(CommandLineInputSpec):
    
    # Inputs
    param_file = File(desc='Parameter .json file', position=-1, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=-1, 
        argstr='', mandatory=True, xor=['param_file'])

    t1_map = File(exists=True, argstr='%s', mandatory=True,
        position=-3, desc='Path to T1 map')
    in_file = File(exists=True, argstr='%s', mandatory=True,
        position=-2, desc='Path to SSFP data')

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
    
class QiDespot2OutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QiDespot2(CommandLine):
    """
    Run DESPOT2 analysis with qidespot2

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot2
    >>> d1 = QiDespot2(prefix='nipype_', param_file='ssfp_params.json')
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot2
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]} }
    >>> d2 = QiDespot2(prefix='nipype_', param_dict=params)
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    """

    _cmd = 'qidespot2'
    input_spec = QiDespot2InputSpec
    output_spec = QiDespot2OutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiDespot2, self)._format_arg(name, spec, value)

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