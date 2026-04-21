#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 12:50:55 2022

@author: avinash
"""
#!/usr/bin/env python

##########
# pptoas #
##########

# pptoas is a command-line program used to simultaneously fit for phases
#    (phis/TOAs), dispersion measures (DMs), frequency**-4 delays (GMs),
#    scattering timescales (taus), and scattering indices (alphas).  Mean flux
#    densities can also be estimated from the fitted model.  Full-functionality
#    is obtained when using pptoas within an interactive python environment.

# Written by Timothy T. Pennucci (TTP; tim.pennucci@nanograv.org).
# Contributions by Scott M. Ransom (SMR) and Paul B. Demorest (PBD).

from __future__ import division
from __future__ import print_function

from builtins import map
from builtins import object
from builtins import range
from builtins import zip

from past.utils import old_div
from types import SimpleNamespace

from mpplib_b35 import *


# cfitsio defines a maximum number of files (NMAXFILES) that can be opened in
# the header file fitsio2.h.  Without calling unload() with PSRCHIVE, which
# touches the archive, I am not sure how to close the files.  So, to avoid the
# loop crashing, set a maximum number of archives for pptoas.  Modern machines
# should be able to handle almost 1000.
max_nfile = 999

# See F0_fact in pplib.py
if F0_fact:
    rm_baseline = True
else:
    rm_baseline = True


class TOA(object):
    """
    TOA class bundles common TOA attributes together with useful functions.
    """

    def __init__(self, archive_band3, archive_band5, frequency, MJD, TOA_error, telescope_band3,
                 telescope_band5,telescope_code_band3,
                 telescope_code_band5, DM=None, DM_error=None, flags={}):
        """
        Form a TOA.

        archive is the string name of the TOA's archive.
        frequency is the reference frequency [MHz] of the TOA.
        MJD is a PSRCHIVE MJD object (the TOA, topocentric).
        TOA_error is the TOA uncertainty [us].
        telescope is the name of the observatory.
        telescope_code is the string written on the TOA line.
        DM is the full DM [cm**-3 pc] associated with the TOA.
        DM_error is the DM uncertainty [cm**-3 pc].
        flags is a dictionary of arbitrary TOA flags
            (e.g. {'subint':0, 'be':'GUPPI'}).
        """
        self.archive_band3 = archive_band3
        self.archive_band5 = archive_band5
        self.frequency = frequency
        self.MJD = MJD
        self.TOA_error = TOA_error
        self.telescope_band3 = telescope_band3
        self.telescope_band5 = telescope_band5
        self.telescope_code_band3 = telescope_code_band3
        self.telescope_code_band5 = telescope_code_band5
        self.DM = DM
        self.DM_error = DM_error
        self.flags = flags
        # This isn't necessary since other functions can just access the dictionary directly
        # for flag in list(flags.keys()):
        #    exec('self.%s = flags["%s"]' % (flag, flag))

    def write_TOA(self, inf_is_zero=True, outfile=None):
        """
        Print a loosely IPTA-formatted TOA to standard output or to file.
        inf_is_zero=True follows the TEMPO/2 convention of writing 0.0 MHz as
            the frequency for infinite-frequency TOAs.
        outfile is the output file name; if None, will print to standard
            output.
        """
        write_TOAs(self, inf_is_zero=inf_is_zero, outfile=outfile, append=True)


class GetTOAs(object):
    """
    GetTOAs is a class with methods to measure TOAs and DMs from data.
    """

    def __init__(self, datafiles_band3, datafiles_band5, modelfile_band3,  modelfile_band5, quiet=False):
        """
        Unpack all of the data and set initial attributes.

        datafiles is either a single PSRCHIVE file name, or a name of a
            metafile containing a list of archive names.
        modelfile is a ppgauss or ppspline model file.  modelfile can also be
            an arbitrary PSRCHIVE archive, although this feature is
            *not*quite*implemented*yet*.
        quiet=True suppresses output.
        """
        if file_is_type(datafiles_band3, "ASCII"):
            self.datafiles_band3 = [datafile_band3[:-1] for datafile_band3 in \
                    open(datafiles_band3, "r").readlines()]
        else:
            self.datafiles_band3 = [datafiles_band3]
        if len(self.datafiles_band3) > max_nfile:
            print("Too many archives.  See/change max_nfile(=%d) in pptoas.py." % max_nfile)
            sys.exit()
            
        if file_is_type(datafiles_band5, "ASCII"):
            self.datafiles_band5 = [datafile_band5[:-1] for datafile_band5 in \
                    open(datafiles_band5, "r").readlines()]
        else:
            self.datafiles_band5 = [datafiles_band5]
        if len(self.datafiles_band5) > max_nfile:
            print ("Too many archives.  See/change max_nfile(=%d) in pptoas.py."%max_nfile)
            sys.exit()
            
        self.is_FITS_model_band3 = file_is_type(modelfile_band3, "FITS")
        self.is_FITS_model_band5 = file_is_type(modelfile_band5, "FITS")
        
        self.modelfile_band3 = modelfile_band3  # the model file in use
        self.modelfile_band5 = modelfile_band5  # the model file in use
# =============================================================================
#         print(self.datafiles_band3)
#         print(self.datafiles_band5)
#         print(self.modelfile_band3)
#         print(self.modelfile_band5)
#         print(ok)
# =============================================================================
        
        self.obs_band3 = []  # observatories from the observations
        self.doppler_fs_band3 = []  # PSRCHIVE Doppler factors from Earth's motion
        self.nu0s_band3 = []  # PSRCHIVE center frequency
        self.obs_band5 = []  # observatories from the observations
        self.doppler_fs_band5 = []  # PSRCHIVE Doppler factors from Earth's motion
        self.nu0s_band5 = []  # PSRCHIVE center frequency
        self.nu_fits = []  # reference frequencies for the fit
        self.nu_refs = []  # reference frequencies for the output
        self.ok_idatafiles_band3 = []  # list of indices for good/examined datafiles
        self.ok_idatafiles_band5 = []  # list of indices for good/examined datafiles
        self.ok_isubs = []  # list of indices for the good subintegrations
        self.epochs_band3 = []  # PSRCHIVE midpoints of the subintegrations
        self.epochs_band5 = []  # PSRCHIVE midpoints of the subintegrations
        self.MJDs = []  # same as epochs, in days
        self.Ps = []  # PSRCHIVE spin period at each epoch
        self.phis = []  # the fitted phase shifts / phi parameter
        self.phi_errs = []  # their uncertainties
        self.TOAs = []  # the fitted TOA
        self.TOA_errs = []  # their uncertainties
        self.DM0s = []  # the stored PSRCHIVE header DMs
        self.DMs = []  # the fitted DMs (may include the Doppler correction)
        self.DM_errs = []  # their uncertainties
        self.DeltaDM_means = []  # fitted single mean DM-DM0
        self.DeltaDM_errs = []  # their uncertainties
        self.GMs = []  # fitted "GM" parameter, from delays that go as nu**-4
        self.GM_errs = []  # their uncertainties
        self.taus = []  # fitted scattering timescales
        self.tau_errs = []  # their uncertainties
        self.alphas = []  # fitted scattering indices
        self.alpha_errs = []  # their uncertainties
        self.scales = []  # fitted per-channel scaling parameters
        self.scale_errs = []  # their uncertainties
        self.snrs = []  # signal-to-noise ratios (S/N)
        self.channel_snrs = []  # per-channel S/Ns
        self.profile_fluxes = []  # estimated per-channel fluxes
        self.profile_flux_errs = []  # their uncertainties
        self.fluxes = []  # estimated overall fluxes
        self.flux_errs = []  # their uncertainties
        self.flux_freqs = []  # their reference frequencies
        self.red_chi2s = []  # reduced chi2 values of the fit
        self.channel_red_chi2s = []  # per-channel reduced chi2 values
        self.covariances = []  # full covariance matrices
        self.nfevals = []  # number of likelihood function evaluations
        self.rcs = []  # return codes from the fit
        self.fit_durations = []  # durations of the fit
        self.order_band3 = []  # order that datafiles are examined (deprecated)
        self.order_band5 = []  # order that datafiles are examined (deprecated)
        self.TOA_list = []  # complete, single list of TOAs
        self.zap_channels = []  # proposed channels to be zapped
        # dictionary of instrumental response characteristics
        self.instrumental_response_dict_band3 = self.ird_band3 = \
            {'DM': 0.0, 'wids': [], 'irf_types': []}
        self.quiet = quiet  # be quiet?

    def get_TOAs(self, datafile_band3=None, datafile_band5=None, tscrunch=False, nu_refs=None, DM0=None,
            bary=True, fit_DM=True, fit_GM=False, fit_scat=False,
            log10_tau=True, scat_guess=None, fix_alpha=False,
            print_phase=True, print_flux=False, print_parangle=False,
            add_instrumental_response=False, addtnl_toa_flags={},
            method='trust-ncg', bounds=None, nu_fits=None, show_plot=False,
            bayes=False, quiet=None, outdir=None):
        """
        Measure TOAs from wideband data accounting for numerous ISM effects.

        datafile defaults to self.datafiles, otherwise it is a single
            PSRCHIVE archive name
        tscrunch=True tscrunches archive before fitting (i.e. make one set of
            measurements per archive)
        nu_refs is a tuple containing two output reference frequencies [MHz],
            one for the TOAs, and the other for the scattering timescales;
            defaults to the zero-covariance frequency between the TOA and DM,
            and the scattering timescale and index, respectively.
        DM0 is the baseline dispersion measure [cm**-3 pc]; defaults to what is
            stored in each datafile.
        bary=True corrects the measured DMs, GMs, taus, and nu_ref_taus based
            on the Doppler motion of the observatory with respect to the solar
            system barycenter.
        fit_DM=False will not fit for DM; if this is the case, you might want
            to set bary to False.
        fit_GM=True will fit for a parameter ('GM') characterizing a delay term
            for each TOA that scales as nu**-4.  Will be highly covariant with
            DM.
        fit_scat=True will fit the scattering timescale and index for each TOA.
        log10_tau=True does the scattering fit with log10(scattering timescale)
            as the parameter.
        scat_guess can be a list of three numbers: a guess of the scattering
            timescale tau [s], its reference frequency [MHz], and a guess of
            the scattering index alpha.  Will be used for all archives;
            supercedes other initial values.
        fix_alpha=True will hold the scattering index fixed, in the case that
            fit_scat==True.  alpha is fixed to the value specified in the
            .gmodel file, or scattering_alpha in pplib.py if no .gmodel is
            provided.
        print_phase=True will print the fitted parameter phi and its
            uncertainty on the TOA line with the flags -phs and -phs_err.
        print_flux=True will print an estimate of the overall flux density and
            its uncertainty on the TOA line.
        print_parangle=True will print the parallactic angle on the TOA line.
        add_instrumental_response=True will account for the instrumental
            response according to the dictionary instrumental_response_dict.
        addtnl_toa_flags are pairs making up TOA flags to be written uniformly
            to all IPTA-formatted TOAs. e.g. ('pta','NANOGrav','version',0.1)
        method is the scipy.optimize.minimize method; currently can be 'TNC',
            'Newton-CG', or 'trust-cng', which are all Newton
            Conjugate-Gradient algorithms.
        bounds is a list of five 2-tuples, giving the lower and upper bounds on
            the phase, dispersion measure, GM, tau, and alpha parameters,
            respectively.  NB: this is only used if method=='TNC'.
        nu_fits is a tuple, analogous to nu_ref, where these reference
            frequencies [MHz] are used in the fit; defaults to a guess at the
            zero-covariance frequency based on signal-to-noise ratios.
        show_plot=True will show a plot of the fitted model, data, and
            residuals at the end of the fitting.  If set to "save" will
            save plots to files.
        quiet=True suppresses output.
        """
        if quiet is None: quiet = self.quiet
        already_warned = False
        warning_message = \
            "You are using an experimental functionality of pptoas!"
        self.nfit = 1
        if fit_DM: self.nfit += 1
        if fit_GM: self.nfit += 1
        if fit_scat: self.nfit += 2
        if fix_alpha: self.nfit -= 1
        self.fit_phi = True
        self.fit_DM = fit_DM
        self.fit_GM = fit_GM
        self.fit_tau = self.fit_alpha = fit_scat
        if fit_scat: self.fit_alpha = not fix_alpha
        self.fit_flags = [int(self.fit_phi), int(self.fit_DM),
                          int(self.fit_GM), int(self.fit_tau), int(self.fit_alpha)]
        self.log10_tau = log10_tau
        if not fit_scat:
            self.log10_tau = log10_tau = False
        if self.fit_GM or fit_scat or self.fit_tau or self.fit_alpha:
            print(warning_message)
            already_warned = True
        self.scat_guess = scat_guess
        nu_ref_tuple = nu_refs
        nu_fit_tuple = nu_fits
        self.DM0 = DM0
        self.bary = bary
        start = time.time()
        tot_duration = 0.0
        
        if datafile_band3 is None:
            datafiles_band3 = self.datafiles_band3
        else:
            datafiles_band3 = [datafile_band3]
        if datafile_band5 is None:
            datafiles_band5 = self.datafiles_band5
        else:
            datafiles_band5 = [datafile_band5]
            
        self.tscrunch = tscrunch
        self.add_instrumental_response = add_instrumental_response
        # print("yes")
        # print(enumerate(datafiles_band3))
        # print(datafiles_band5)
# =============================================================================
#         print("avinash: datafiles_band3",datafiles_band3)
#         print("avinash: datafiles_band5",datafiles_band5)
#         print(ok)
# =============================================================================
        
        for iarch in range(len(datafiles_band3)):
            datafile_band3 = datafiles_band3[iarch]
            datafile_band5 = datafiles_band5[iarch]
# =============================================================================
#             print(iarch)
#             print(datafile_band3)
#             print(datafile_band5)
#             print(ok)
# =============================================================================
            fit_duration = 0.0
            # Load data
            try:
                data_band3 = load_data(datafile_band3, dedisperse=False,
                                 dededisperse=False, tscrunch=tscrunch,
                                 pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                                 flux_prof=False, refresh_arch=False, return_arch=False,
                                 quiet=quiet)
                if data_band3.dmc:
                    if not quiet:
                        print("%s is dedispersed (dmc = 1).  Reloading it." % \
                              datafile_band3)
                    # continue
                    data_band3 = load_data(datafile_band3, dedisperse=False,
                                     dededisperse=True, tscrunch=tscrunch,
                                     pscrunch=True, fscrunch=False,
                                     rm_baseline=rm_baseline, flux_prof=False,
                                     refresh_arch=False, return_arch=False, quiet=quiet)
                if not len(data_band3.ok_isubs):
                    if not quiet:
                        print("No subints to fit for %s.  Skipping it." % \
                              datafile_band3)
                    continue
                else:
                    self.ok_idatafiles_band3.append(iarch)
            except RuntimeError:
                if not quiet:
                    print("Cannot load_data(%s).  Skipping it." % datafile_band3)
                continue
            
            try:
                data_band5 = load_data(datafile_band5, dedisperse=False,
                                 dededisperse=False, tscrunch=tscrunch,
                                 pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                                 flux_prof=False, refresh_arch=False, return_arch=False,
                                 quiet=quiet)
                if data_band5.dmc:
                    if not quiet:
                        print("%s is dedispersed (dmc = 1).  Reloading it." % \
                              datafile_band5)
                    # continue
                    data_band5 = load_data(datafile_band5, dedisperse=False,
                                     dededisperse=True, tscrunch=tscrunch,
                                     pscrunch=True, fscrunch=False,
                                     rm_baseline=rm_baseline, flux_prof=False,
                                     refresh_arch=False, return_arch=False, quiet=quiet)
                if not len(data_band5.ok_isubs):
                    if not quiet:
                        print("No subints to fit for %s.  Skipping it." % \
                              datafile_band5)
                    continue
                else:
                    self.ok_idatafiles_band5.append(iarch)
            except RuntimeError:
                if not quiet:
                    print("Cannot load_data(%s).  Skipping it." % datafile_band5)
                continue
            # print("yes")
            # Unpack the data dictionary into the local namespace; see load_data
            # for dictionary keys.
            # for key in list(data.keys()):
            #     exec(key + " = data['" + key + "']")
            # BWM: Python 3 no longer supports updating local namespace via exec statement/function
            # Also, since we're setting local variables, it's not safe to just update the object __dict__
            # So the best option is to instead load the dictionary into a SimpleNamespace and
            # access key/values as you would an object attribute
            
# =============================================================================
#             print(data_band3)
#             print(data_band5)
#             print(ok)
# =============================================================================
# =============================================================================
# # =============================================================================
# #             print(data_band3)
# #             print(data_band5)
# #             print(ok)
# # =============================================================================
# =============================================================================
            d3 = SimpleNamespace(**data_band3)
            
            nsub_band3 = d3.nsub
            nchan_band3 = d3.nchan
            nbin_band3 = d3.nbin
            npol_band3 = d3.npol
            d5 = SimpleNamespace(**data_band5)
            nsub_band5 = d5.nsub
            nchan_band5 = d5.nchan
            nbin_band5 = d5.nbin
            npol_band5 = d5.npol
            
# =============================================================================
#             print(nsub_band3,nchan_band3,nbin_band3,npol_band3)
#             print(nsub_band5,nchan_band5,nbin_band5,npol_band5)
#             print(ok)
# =============================================================================
            if d3.source is None: d3.source = "noname"
            if d5.source is None: d5.source = "noname"
            # Observation info
            obs_band3 = DataBunch(telescope_band3=d3.telescope, backend_band3=d3.backend,
                            frontend_band3=d3.frontend)
            nu_fits = list(np.zeros([(nsub_band3), 3], dtype=np.float64))
            nu_refs = list(np.zeros([(nsub_band3), 3], dtype=np.float64))
            phis = np.zeros(nsub_band3, dtype=np.double)
            phi_errs = np.zeros(nsub_band3, dtype=np.double)
            TOAs = np.zeros(nsub_band3, dtype="object")
            TOA_errs = np.zeros(nsub_band3, dtype="object")
            DMs = np.zeros(nsub_band3, dtype=np.float64)
            DM_errs = np.zeros(nsub_band3, dtype=np.float64)
            GMs = np.zeros(nsub_band3, dtype=np.float64)
            GM_errs = np.zeros(nsub_band3, dtype=np.float64)
            taus = np.zeros(nsub_band3, dtype=np.float64)
            tau_errs = np.zeros(nsub_band3, dtype=np.float64)
            alphas = np.zeros(nsub_band3, dtype=np.float64)
            alpha_errs = np.zeros(nsub_band3, dtype=np.float64)
            scales = np.zeros([nsub_band3, nchan_band3+nchan_band5], dtype=np.float64)
            scale_errs = np.zeros([nsub_band3, nchan_band3+nchan_band5], dtype=np.float64)
            snrs = np.zeros(nsub_band3, dtype=np.float64)
            channel_snrs = np.zeros([nsub_band3, nchan_band3+nchan_band5], dtype=np.float64)
            profile_fluxes = np.zeros([nsub_band3, nchan_band3+nchan_band5], dtype=np.float64)
            profile_flux_errs = np.zeros([nsub_band3, nchan_band3+nchan_band5], dtype=np.float64)
            fluxes = np.zeros(nsub_band3, dtype=np.float64)
            flux_errs = np.zeros(nsub_band3, dtype=np.float64)
            flux_freqs = np.zeros(nsub_band3, dtype=np.float64)
            red_chi2s = np.zeros(nsub_band3, dtype=np.float64)
            covariances = np.zeros([nsub_band3, self.nfit, self.nfit],
                                   dtype=np.float64)
            nfevals = np.zeros(nsub_band3, dtype="int")
            rcs = np.zeros(nsub_band3, dtype="int")
            
            obs_band5 = DataBunch(telescope_band5=d5.telescope, backend_band5=d5.backend,
                            frontend_band5=d5.frontend)
            ok_ichans_b35 = np.zeros([nsub_band3,nchan_band3+nchan_band5], dtype="int")
            
            # print("ok")
            # PSRCHIVE epochs are *midpoint* of the integration
            MJDs_band3 = np.array([d3.epochs[isub].in_days() \
                             for isub in range(nsub_band3)], dtype=np.double)
            DM_stored_band3 = d3.DM  # same as = arch.get_dispersion_measure()
            MJDs_band5 = np.array([d5.epochs[isub].in_days() \
                             for isub in range(nsub_band5)], dtype=np.double)
            DM_stored_band5 = d5.DM  # same as = arch.get_dispersion_measure()
# =============================================================================
#             print("avinash: MJD_band3",MJDs_band3)
#             print("avinash: MJD_band5",MJDs_band5)
#             print("avinash: DM_band3",DM_stored_band3)
#             print("avinash: DM_band5",DM_stored_band5)
#             print(ok)
# =============================================================================
            if self.DM0 is None:
                DM0 = DM_stored_band3
            else:
                DM0 = self.DM0
            if self.is_FITS_model_band3:
                print("if self.is_FITS_model_band3: is not executed for band 3 model as it is false")
# =============================================================================
#                 if not already_warned:
#                     print(warning_message)
#                     already_warned = True
#                 model_data_band3 = load_data(self.modelfile_band3, dedisperse=False,
#                     dededisperse=False, tscrunch=True, pscrunch=True,
#                     fscrunch=False, rm_baseline=True, flux_prof=False,
#                     #fscrunch=False, rm_baseline=False, flux_prof=False,
#                     refresh_arch=False, return_arch=False, quiet=True)
#                 model_band3 = (model_data_band3.masks * model_data_band3.subints)[0,0]
#                 if model_data_band3.nbin != nbin:
#                     print("Model_band3 nbin %d != data nbin %d for %s; skipping it."\
#                             %(model_data_band3.nbin, nbin_band3, datafile_band3))
#                     continue
#                 if model_data_band3.nchan == 1:
#                     model_band3 = np.tile(model_band3[0], len(freqs[0])).reshape(
#                             len(freqs[0]), nbin)
#                 if model_data_band3.nchan != nchan_band3:
#                     print("Model_band3 nchan %d != data nchan %d for %s; skipping it."%(model_data_band3.nchan, nchan_band3, datafile_band3))
#                     continue
# =============================================================================
            if self.is_FITS_model_band5:
                print("if self.is_FITS_model_band3: is not executed for band 3 model as it is false")
# =============================================================================
#                 if not already_warned:
#                     print(warning_message)
#                     already_warned = True
#                 model_data_band5 = load_data(self.modelfile_band5, dedisperse=False,
#                     dededisperse=False, tscrunch=True, pscrunch=True,
#                     fscrunch=False, rm_baseline=True, flux_prof=False,
#                     #fscrunch=False, rm_baseline=False, flux_prof=False,
#                     refresh_arch=False, return_arch=False, quiet=True)
#                 model_band5 = (model_data_band5.masks * model_data_band5.subints)[0,0]
#                 if model_data_band5.nbin != nbin:
#                     print("Model_band5 nbin %d != data nbin %d for %s; skipping it."\
#                             %(model_data_band5.nbin, nbin_band5, datafile_band5))
#                     continue
#                 if model_data_band5.nchan == 1:
#                     model_band5 = np.tile(model_band5[0], len(freqs[0])).reshape(
#                             len(freqs[0]), nbin)
#                 if model_data_band5.nchan != nchan_band5:
#                     print("Model_band5 nchan %d != data nchan %d for %s; skipping it."%(model_data_band5.nchan, nchan_band5, datafile_band5))
#                     continue
# =============================================================================

                
            if not quiet:
                print("\nEach of the %d TOAs is approximately %.2f s" % (
                        len(d3.ok_isubs), old_div(d3.integration_length, nsub_band3)))
                print("\nEach of the %d TOAs is approximately %.2f s" % (
                    len(d5.ok_isubs), old_div(d5.integration_length, nsub_band5)))
                
            itoa = 1
                        
            
            for isub in d3.ok_isubs:
                if self.is_FITS_model_band3 and \
                        np.any(model_data_band3.freqs[0] - d3.freqs[isub]): # tscrunched
                            print("Frequency mismatch between template and data in band3!")
                sub_id_band3 = datafile_band3 + "_%d"%isub
                epoch_band3 = d3.epochs[isub]
                MJD_band3 = MJDs_band3[isub]
                P_band3 = d3.Ps[isub]
# =============================================================================
#                 print("avinash: sub_id_band3",sub_id_band3)
#                 print("avinash: epoch_band3",epoch_band3)
#                 print("avinash: MJD_band3",MJD_band3)
#                 print("avinash: P_band3",P_band3)
#                 print(ok)
# =============================================================================
                
                if not self.is_FITS_model_band3:
                    # Read model
                    try:
                        if not fit_scat:
                            self.model_name_band3, self.ngauss_band3, model_band3 = read_model(
                                self.modelfile_band3, d3.phases, d3.freqs[isub],
                                d3.Ps[isub],
                                quiet=bool(quiet + (itoa - 1)))
        #                 else:
                #                 self.model_name_band3, self.ngauss_band3, full_model_band3 = \
                #                     read_model(self.modelfile_band3, d3.phases,
                #                                d3.freqs[isub], d3.Ps[isub],
                #                                quiet=bool(quiet + (itoa - 1)))
                #                 self.model_name_band3, self.model_code_band3, \
                #                 self.model_nu_ref, self.ngauss, \
                #                 self.gparams, model_fit_flags, self.alpha, \
                #                 model_fit_alpha = read_model(
                #                     self.modelfile_band3,
                #                     quiet=bool(quiet + (itoa - 1)))
                #                 unscat_params = np.copy(self.gparams)
                #                 unscat_params[1] = 0.0
                #                 model = unscat_model = gen_gaussian_portrait(
                #                     self.model_code, unscat_params, 0.0,
                #                     d.phases, d.freqs[isub], self.model_nu_ref)
                    except (UnboundLocalError,UnicodeDecodeError):
                        self.model_name_band3, model_band3 = read_spline_model(
                            self.modelfile_band3, d3.freqs[isub], nbin_band3,
                            quiet=True)  # bool(quiet+(itoa-1)))
                # else:
                ##THESE FREQUENCIES WILL BE OFF IF AVERAGED CHANNELS##
                #    print model_data.freqs[0, ok_ichans[isub]] - \
                #            freqs[isub,ok_ichans[isub]]
                freqsx_band3 = d3.freqs[isub, d3.ok_ichans[isub]]
                weightsx_band3 = d3.weights[isub, d3.ok_ichans[isub]]
                portx_band3 = d3.subints[isub, 0, d3.ok_ichans[isub]]
                modelx_band3 = model_band3[d3.ok_ichans[isub]]
                if add_instrumental_response and \
                        (self.ird_band3['DM'] or len(self.ird_band3['wids'])):
                    print("band3 instrumental response is is not taken care of please edit code for this")
# =============================================================================
#                     inst_resp_port_FT_band3 = instrumental_response_port_FT(
#                         nbin_band3, freqsx_band3, self.ird_band3['DM'], P_band3,
#                         self.ird_band3['wids'], self.ird_band3['irf_types'])
#                     modelx_band3 = fft.irfft(inst_resp_port_FT_band3 * \
#                                        fft.rfft(modelx_band3, axis=-1), axis=-1)
# =============================================================================
                SNRsx_band3 = d3.SNRs[isub, 0, d3.ok_ichans[isub]]
                # NB: Time-domain uncertainties below
                errs_band3 = d3.noise_stds[isub, 0, d3.ok_ichans[isub]]
                # nu_fit is the reference frequency for parameters in the fit
                nu_mean_band3 = freqsx_band3.mean()
                
                if self.is_FITS_model_band5 and \
                        np.any(model_data_band5.freqs[0] - d5.freqs[isub]): # tscrunched
                            print("Frequency mismatch between template and data in band5!")
                sub_id_band5 = datafile_band5 + "_%d"%isub
                epoch_band5 = d5.epochs[isub]
                MJD_band5 = MJDs_band5[isub]
                P_band5 = d5.Ps[isub]
# =============================================================================
#                 print("avinash: sub_id_band5",sub_id_band5)
#                 print("avinash: epoch_band5",epoch_band5)
#                 print("avinash: MJD_band5",MJD_band5)
#                 print("avinash: P_band5",P_band5)
#                 print(ok)
# =============================================================================
                if not self.is_FITS_model_band5:
                    # Read model
                    try:
                        if not fit_scat:
                            self.model_name_band5, self.ngauss_band5, model_band5 = read_model(
                                self.modelfile_band5, d5.phases, d5.freqs[isub],
                                d5.Ps[isub],
                                quiet=bool(quiet + (itoa - 1)))
        #                 else:
                #                 self.model_name_band5, self.ngauss_band5, full_model_band5 = \
                #                     read_model(self.modelfile_band5, d5.phases,
                #                                d5.freqs[isub], d5.Ps[isub],
                #                                quiet=bool(quiet + (itoa - 1)))
                #                 self.model_name_band5, self.model_code_band5, \
                #                 self.model_nu_ref, self.ngauss, \
                #                 self.gparams, model_fit_flags, self.alpha, \
                #                 model_fit_alpha = read_model(
                #                     self.modelfile_band5,
                #                     quiet=bool(quiet + (itoa - 1)))
                #                 unscat_params = np.copy(self.gparams)
                #                 unscat_params[1] = 0.0
                #                 model = unscat_model = gen_gaussian_portrait(
                #                     self.model_code, unscat_params, 0.0,
                #                     d.phases, d.freqs[isub], self.model_nu_ref)
                    except (UnboundLocalError,UnicodeDecodeError):
                        self.model_name_band5, model_band5 = read_spline_model(
                            self.modelfile_band5, d5.freqs[isub], nbin_band5,
                            quiet=True)  # bool(quiet+(itoa-1)))
                # else:
                ##THESE FREQUENCIES WILL BE OFF IF AVERAGED CHANNELS##
                #    print model_data.freqs[0, ok_ichans[isub]] - \
                #            freqs[isub,ok_ichans[isub]]
                freqsx_band5 = d5.freqs[isub, d5.ok_ichans[isub]]
                weightsx_band5 = d5.weights[isub, d5.ok_ichans[isub]]
                portx_band5 = d5.subints[isub, 0, d5.ok_ichans[isub]]
                modelx_band5 = model_band5[d5.ok_ichans[isub]]
                if add_instrumental_response and \
                        (self.ird_band5['DM'] or len(self.ird_band5['wids'])):
                    print("band3 instrumental response is is not taken care of please edit code for this")
# =============================================================================
#                     inst_resp_port_FT_band5 = instrumental_response_port_FT(
#                         nbin_band5, freqsx_band5, self.ird_band5['DM'], P_band5,
#                         self.ird_band5['wids'], self.ird_band5['irf_types'])
#                     modelx_band5 = fft.irfft(inst_resp_port_FT_band5 * \
#                                        fft.rfft(modelx_band5, axis=-1), axis=-1)
# =============================================================================
                SNRsx_band5 = d5.SNRs[isub, 0, d5.ok_ichans[isub]]
                # NB: Time-domain uncertainties below
                errs_band5 = d5.noise_stds[isub, 0, d5.ok_ichans[isub]]
                # nu_fit is the reference frequency for parameters in the fit
                nu_mean_band5 = freqsx_band5.mean()
                

                ##Needs to check for J1643 or something which has profile evolution
# =============================================================================
#                 print("avinash: self.model_name_band3",self.model_name_band3)
#                 print("avinash: model_shape_band3",np.shape(model_band3))
#                 print("avinash: model_band3",model_band3)
#                 print("avinash: self.model_name_band5",self.model_name_band5)
#                 print("avinash: model_shape_band5",np.shape(model_band5))
#                 print("avinash: model_band5",model_band5)
#                 print(ok)
# =============================================================================
# =============================================================================
#                 print("avinash: freqsx_band3",freqsx_band3)
#                 print("avinash: weightsx_band3",weightsx_band3)
#                 print("avinash: portx_band3",portx_band3)
#                 print("avinash: modelx_band3",modelx_band3)                
#                 print("avinash: freqsx_band5",freqsx_band5)
#                 print("avinash: weightsx_band5",weightsx_band5)
#                 print("avinash: portx_band5",portx_band5)
#                 print("avinash: modelx_band5",modelx_band5)
#                 print(ok)
# =============================================================================

# =============================================================================
#                 print("avinash: SNRsx_band3",SNRsx_band3)
#                 print("avinash: errs_band3",errs_band3)
#                 print("avinash: nu_mean_band3",nu_mean_band3)
#                 print("avinash: SNRsx_band5",SNRsx_band5)
#                 print("avinash: errs_band5",errs_band5)
#                 print("avinash: nu_mean_band5",nu_mean_band5)
#                 print(ok)
# =============================================================================
                
                if nu_fit_tuple is None:
                    # NB: the subints are dedispersed at different nu_fit.
                    freqsx = np.concatenate((freqsx_band3,freqsx_band5))
                    SNRsx = np.concatenate((SNRsx_band3,SNRsx_band5))
                    
                    nu_fit = guess_fit_freq(freqsx, SNRsx)
                    nu_fit_DM = nu_fit_GM = nu_fit_tau = nu_fit
# =============================================================================
#                     print("avinash: freqsx",freqsx)
#                     print("avinash: SNRsx",SNRsx)
#                     print("avinash: nu_fit",nu_fit)
#                     print(ok)
# =============================================================================
                else:
                    nu_fit_DM = nu_fit_GM = nu_fit_tuple[0]
                    nu_fit_tau = nu_fit_tuple[-1]
                print(nu_fit_DM, nu_fit_GM, nu_fit_tau)
                nu_fits[isub] = [nu_fit_DM, nu_fit_GM, nu_fit_tau]
                if nu_ref_tuple is None:
                    nu_ref = None
                    nu_ref_DM = nu_ref_GM = nu_ref_tau = nu_ref
                else:
                    nu_ref_DM = nu_ref_GM = nu_ref_tuple[0]
                    nu_ref_tau = nu_ref_tuple[-1]
                    if bary and nu_ref_tau:  # from bary to topo below
                        print("if bary and nu_ref_tau needs to be updated. Currently unavailable due to lack of knowledge on implementing to dual bands")
        #                 nu_ref_tau /= d.doppler_factors[isub]
                nu_refs[isub] = [nu_ref_DM, nu_ref_GM, nu_ref_tau]
                
# =============================================================================
#                 print("avinash: nu_refs[isub]",nu_refs[isub])
#                 print(ok)
# =============================================================================
                

                ###################
                # INITIAL GUESSES #
                ###################
                nu_mean = freqsx.mean()
                DM_guess = DM_stored_band3
                
                #print("avinash: epoch_band3",epoch_band3)
                #print("avinash: epoch_band5",epoch_band5)
                #print("avinash: MJD_band3",MJD_band3)
                #print("avinash: MJD_band5",MJD_band5)
                #print("avinash: P_band3",P_band3)
                #print("avinash: P_band5",P_band5)
                #print("avinash: d3.backend_delay",d3.backend_delay)
                #print("avinash: d5.backend_delay",d5.backend_delay)
# =============================================================================
#                 results.TOA = epoch_band3 + pr.MJD(
#                     old_div(((results.phi * P_band3) + d3.backend_delay),
#                             (3600 * 24.)))
# =============================================================================
                
                phase_offset = old_div(((MJD_band3 - MJD_band5)*24.*3600 + (d3.backend_delay - d5.backend_delay)),P_band5)
                
                print("phase_offset",phase_offset)
                
                
                rot_port_band3 = rotate_data(portx_band3, 0.0, DM_guess, P_band3, freqsx_band3,
                                       nu_mean_band3)  # why not nu_fit?
                rot_port_band5 = rotate_data(portx_band5, phase_offset, DM_guess, P_band5, freqsx_band5,
                                       nu_mean_band3)  # why not nu_fit?
                rot_prof_band3 = np.average(rot_port_band3, axis=0, weights=weightsx_band3)
                rot_prof_band5 = np.average(rot_port_band5, axis=0, weights=weightsx_band5)
                GM_guess = 0.0
                tau_guess = 0.0
                alpha_guess = 0.0
                
# =============================================================================
#                 print("avinash: portx_band3",portx_band3)
#                 print("avinash: portx_band5",portx_band5)
# =============================================================================
# =============================================================================
#                 print("avinash: DM_guess",DM_guess)
#                 print("avinash: rot_port_band3",rot_port_band3)
#                 print("avinash: rot_port_band5",rot_port_band5)
#                 print("avinash: rot_prof_band3",rot_prof_band3)
#                 print("avinash: rot_prof_band5",rot_prof_band5)
# =============================================================================
                
#=============================================================================
# =============================================================================
#                 plt.plot(figsize=(20,20))
#                 plt.imshow(portx_band3,extent=(0,1,0,32), aspect='auto')
#                 plt.title("portx_34")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
#                 plt.plot(figsize=(20,20))
#                 plt.imshow(portx_band5,extent=(0,1,0,32), aspect='auto')
#                 plt.title("portx_45")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
#                 plt.plot(figsize=(20,20))
#                 plt.imshow(modelx_band3,extent=(0,1,0,32), aspect='auto')
#                 plt.title("modelx_34")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
#                 plt.plot(figsize=(20,20))
#                 plt.imshow(modelx_band5,extent=(0,1,0,32), aspect='auto')
#                 plt.title("modelx_45")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
#                 plt.plot(figsize=(20,20))
#                 plt.imshow(rot_port_band3,extent=(0,1,0,64), aspect='auto')
#                 plt.title("rot_port_34")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
#                 plt.plot(figsize=(20,20))
#                 plt.imshow(rot_port_band5,extent=(0,1,0,64), aspect='auto')
#                 plt.title("rot_port_45")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
#                 
#                 print(np.where(rot_prof_band3==max(rot_prof_band3)))
#                 print(np.where(rot_prof_band5==max(rot_prof_band5)))
#                 plt.plot(figsize=(20,20))
#                 plt.plot(np.arange(len(rot_prof_band3)),rot_prof_band3)
#                 plt.title("rot_prof_34")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
#                 plt.plot(figsize=(20,20))
#                 plt.plot(np.arange(len(rot_prof_band5)),rot_prof_band5)
#                 plt.title("rot_prof_45")
#                 plt.xlabel('nbin')
#                 plt.ylabel('nchan')
#                 plt.show()
# # =============================================================================
#                 print(ok)
# =============================================================================
# =============================================================================
#=============================================================================
                
                if fit_scat:
                    print("fit_scat is not available in this dual band version")
# =============================================================================
#                     if self.scat_guess is not None:
#                         tau_guess_s, tau_guess_ref, alpha_guess = self.scat_guess
#                         tau_guess = (old_div(tau_guess_s, P)) * \
#                                     (old_div(nu_fit_tau, tau_guess_ref)) ** alpha_guess
#                     else:
#                         if hasattr(self, 'alpha'):
#                             alpha_guess = self.alpha
#                         else:
#                             alpha_guess = scattering_alpha
#                         if hasattr(self, 'gparams'):
#                             tau_guess = (old_div(self.gparams[1], P)) * \
#                                         (old_div(nu_fit_tau, self.model_nu_ref)) ** alpha_guess
#                         else:
#                             tau_guess = 0.0  # nbin**-1?
#                             # tau_guess = guess_tau(...)
#                     model_prof_scat = fft.irfft(scattering_portrait_FT(
#                         np.array([scattering_times(tau_guess, alpha_guess,
#                                                    nu_fit_tau, nu_fit_tau)]), nbin)[0] * fft.rfft(
#                         modelx.mean(axis=0)))
#                     phi_guess = fit_phase_shift(rot_prof,
#                                                 model_prof_scat, Ns=100).phase
#                     if self.log10_tau:
#                         if tau_guess == 0.0: tau_guess = nbin ** -1
#                         tau_guess = np.log10(tau_guess)
# =============================================================================
                else:
                    # NB: Ns should be larger than nbin for very low S/N data,
                    # especially in the case of noisy models...
                    
                    phi_guess_band3 = fit_phase_shift(rot_prof_band3, modelx_band3.mean(axis=0),
                                                Ns=100).phase
                    phi_guess_band5 = fit_phase_shift(rot_prof_band5, modelx_band5.mean(axis=0),
                                                Ns=100).phase
# =============================================================================
#                     print("avinash: phi_guess_band3",phi_guess_band3)
#                     print("avinash: phi_guess_band5",phi_guess_band5)
#                     print("avinash: nu_mean_band3",nu_mean_band3)
#                     print("aviansh: nu_mean_band5",nu_mean_band5)
#                     print("avinash: nu_fit_DM", nu_fit_DM)
# =============================================================================
                    
                phi_guess = phase_transform(phi_guess_band3, DM_guess, nu_mean_band3,
                                            nu_fit_DM, P_band3, mod=True)  # why not use nu_fit at first?
                print(phi_guess)
                # Need a status bar?
                param_guesses = [phi_guess, DM_guess, GM_guess, tau_guess,
                                 alpha_guess]
# =============================================================================
#                 print("avinash: param_guess",param_guesses)
# =============================================================================
                
                if bounds is None and method == 'TNC':
                    phi_bounds = (None, None)
                    DM_bounds = (None, None)
                    GM_bounds = (None, None)
                    if not self.log10_tau:
                        tau_bounds = (0.0, None)
                    else:
                        tau_bounds = (np.log10((10 * nbin) ** -1), None)
                    alpha_bounds = (-10.0, 10.0)
                    bounds = [phi_bounds, DM_bounds, GM_bounds, tau_bounds,
                              alpha_bounds]
                
                ###########
                # THE FIT #
                ###########
                if not quiet: print("Fitting for TOA #%d" % (itoa))
                if len(freqsx) == 1:
                    ##For combination this thing needs to be edited in detail
                    fit_flags = [1, 0, 0, 0, 0]
                    if not quiet:
                        print("TOA #%d only has 1 frequency channel...fitting for phase only..." % (itoa))
                elif len(freqsx) == 2 and self.fit_DM and self.fit_GM:
                    # prioritize DM fit
                    fit_flags[2] = 0
                    if not quiet:
                        print("TOA #%d only has 2 frequency channels...fitting for phase and DM only..." % (itoa))
                else:
                    fit_flags = list(np.copy(self.fit_flags))
# =============================================================================
#                 print("avinash: nu_ref", nu_refs[isub])
#                 print(ok)
# =============================================================================
#               avinash: Additional phase offset added as parameter
################################################################################
                filename = datafile_band3
###########################################################################
                if not mlan:
                    results = fit_portrait_full(portx_band3, portx_band5, modelx_band3, modelx_band5,
                                                param_guesses, P_band3, phase_offset,
                                                freqsx_band3, freqsx_band5, nu_fits[isub], nu_refs[isub], 
                                                errs_band3, errs_band5, fit_flags,
                                                bounds, self.log10_tau, option=0, 
                                                sub_id_band3=sub_id_band3, sub_id_band5=sub_id_band5, 
                                                method=method, is_toa=True, bayes=bayes, quiet=quiet, 
                                                filename = filename, pool = pool, outdir=outdir)
                if mlan:
                    results = fit_portrait_full_mlan(portx_band3, portx_band5, modelx_band3, modelx_band5,
                                                param_guesses, P_band3, phase_offset,
                                                freqsx_band3, freqsx_band5, nu_fits[isub], nu_refs[isub], 
                                                errs_band3, errs_band5, fit_flags,
                                                bounds, self.log10_tau, option=0, 
                                                sub_id_band3=sub_id_band3, sub_id_band5=sub_id_band5, 
                                                method=method, is_toa=True, bayes=bayes, quiet=quiet, 
                                                filename = filename, pool = pool, outdir=outdir)
                # Old code
                # results = fit_portrait(portx, modelx,
                #        np.array([phi_guess, DM_guess]), P, freqsx,
                #        nu_fit_DM, nu_ref_DM, errs, bounds=bounds, id=sub_id,
                #        quiet=quiet)
                # results.phi = results.phase
                # results.phi_err = results.phase_err
                # results.GM = results.GM_err = None
                # results.tau = results.tau_err = None
                # results.alpha = results.alpha_err = None
                # results.covariance_matrix = np.zeros([2,2])
                # results.nu_DM = results.nu_GM = results.nu_tau =results.nu_ref
                # Old code for fitting just phase...
                # else:  #1-channel hack
                #    if not quiet:
                #        print "TOA only has %d frequency channel!..."%len(
                #                freqsx)
                #        print "...using Fourier phase gradient routine to fit phase only..."
                #    results = fit_phase_shift(portx[0], modelx[0], errs[0],
                #            Ns=nbin)
                #    results.phi = results.phase
                #    results.phi_err = results.phase_err
                #    results.DM = results.DM_err = None
                #    results.GM = results.GM_err = None
                #    results.tau = results.tau_err = None
                #    results.alpha = results.alpha_err = None
                #    results.nu_DM, results.nu_GM, results.nu_tau = \
                #            [freqsx[0], freqsx[0], freqsx[0]]
                #    results.nfeval = 0
                #    results.return_code = -2
                #    results.scales = np.array([results.scale])
                #    results.scale_errs = np.array([results.scale_err])
                #    results.covariance_matrix = np.identity(self.nfit)
                fit_duration += results.duration

                ####################
                #  CALCULATE  TOA  #
                ####################
            
                results.TOA = epoch_band3 + pr.MJD(
                    old_div(((results.phi * P_band3) + d3.backend_delay),
                            (3600 * 24.)))
                results.TOA_err = results.phi_err * P_band3 * 1e6  # [us]
                
                ############################################################################
                results.TOA_inf = epoch_band3 + pr.MJD(
                    old_div(((results.phi_inf * P_band3) + d3.backend_delay),
                            (3600 * 24.)))
                results.TOA_inf_err = results.phi_inf_err * P_band3 * 1e6  # [us]
                ################################################################################

                ######################
                # DOPPLER CORRECTION #
                ######################
                # This correction should fix Doppler-induced annual variations
                # to DM(t), but will not fix Doppler-induced /orbital/
                # variations to DM(t).

                if self.bary:  # Default is True
                    df = d3.doppler_factors[isub]
                    if fit_flags[1]:
                        # NB: The following eqution was incorrectly reversed in
                        #     the original paper Pennucci et al. (2014),
                        #     printed as DM_bary = DM_topo / df.
                        results.DM *= df  # NB: No longer the *fitted* value!
                    if fit_flags[2]:
                        results.GM *= df ** 3  # NB: No longer the *fitted* value!
                else:
                    df = 1.0
                

                #################
                # ESTIMATE FLUX #
                #################
                
                if print_flux:
                    print("Flux estimentation is not implemented")
# =============================================================================
#                     if results.tau != 0.0:
#                         if self.log10_tau:
#                             tau = 10 ** results.tau
#                         else:
#                             tau = results.tau
#                         alpha = results.alpha
#                         
#                         scat_model = fft.irfft(scattering_portrait_FT(
#                             scattering_times(tau, alpha, freqsx,
#                                              results.nu_tau), data.nbin) * \
#                                                fft.rfft(modelx, axis=1), axis=1)
#                     else:
#                         scat_model = np.copy(modelx)
#                     scat_model_means = scat_model.mean(axis=1)
#                     profile_fluxes[isub, d.ok_ichans[isub]] = scat_model_means * \
#                                                               results.scales
#                     profile_flux_errs[isub, d.ok_ichans[isub]] = abs(
#                         scat_model_means) * results.scale_errs
#                     flux, flux_err = weighted_mean(profile_fluxes[isub,
#                                                                   d.ok_ichans[isub]], profile_flux_errs[isub,
#                                                                                                         d.ok_ichans[
#                                                                                                             isub]])
#                     flux_freq, flux_freq_err = weighted_mean(freqsx,
#                                                              profile_flux_errs[isub, d.ok_ichans[isub]])
#                     fluxes[isub] = flux
#                     flux_errs[isub] = flux_err
#                     flux_freqs[isub] = flux_freq
# 
# =============================================================================
                
                nu_refs[isub] = [results.nu_DM, results.nu_GM, results.nu_tau]
                phis[isub] = results.phi
                phi_errs[isub] = results.phi_err
                TOAs[isub] = results.TOA
                TOA_errs[isub] = results.TOA_err
                DMs[isub] = results.DM
                DM_errs[isub] = results.DM_err
                GMs[isub] = results.GM
                GM_errs[isub] = results.GM_err
                taus[isub] = results.tau
                tau_errs[isub] = results.tau_err
                alphas[isub] = results.alpha
                alpha_errs[isub] = results.alpha_err
                nfevals[isub] = results.nfeval
                rcs[isub] = results.return_code
                ok_ichans_b35[isub] = np.arange(len(d3.ok_ichans[isub])+len(d5.ok_ichans[isub]))
                print(ok_ichans_b35[isub])
                scales[isub, ok_ichans_b35[isub]] = results.scales
                scale_errs[isub, ok_ichans_b35[isub]] = results.scale_errs
                snrs[isub] = results.snr
                channel_snrs[isub, ok_ichans_b35[isub]] = results.channel_snrs
                try:
                    covariances[isub] = results.covariance_matrix
                except ValueError:
                    for ii, ifit in enumerate(np.where(fit_flags)[0]):
                        for jj, jfit in enumerate(np.where(fit_flags)[0]):
                            covariances[isub][ifit, jfit] = \
                                results.covariance_matrix[ii, jj]
                red_chi2s[isub] = results.red_chi2
                # Compile useful TOA flags
                # Add doppler_factor?
                #print(df)
            
                toa_flags = {}
                if not fit_flags[1]:
                    results.DM = None
                    results.DM_err = None
                if fit_flags[2]:
                    toa_flags['gm'] = results.GM
                    toa_flags['gm_err'] = results.GM_err
                if fit_flags[3]:
                    if self.log10_tau:
                        toa_flags['scat_time'] = old_div(10 ** results.tau * P_band3, df) * 1e6
                        # usec, w/ df
                        toa_flags['log10_scat_time'] = results.tau + \
                                                       np.log10(old_div(P_band3, df))  # w/ df
                        toa_flags['log10_scat_time_err'] = results.tau_err
                    else:
                        toa_flags['scat_time'] = old_div(results.tau * P_band3, df) * 1e6
                        # usec, w/ df
                        toa_flags['scat_time_err'] = old_div(results.tau_err * P_band3, df) \
                                                     * 1e6  # usec, w/ df
                    toa_flags['scat_ref_freq'] = results.nu_tau * df  # w/ df
                    toa_flags['scat_ind'] = results.alpha
                if fit_flags[4]:
                    toa_flags['scat_ind_err'] = results.alpha_err

                toa_flags['be'] = d3.backend + "_&_" + d5.backend
                toa_flags['fe'] = d3.frontend + "_&_" + d5.frontend
                toa_flags['f'] = d3.frontend + "_" + d3.backend + "_" + str(int(d3.bw)) + "_&_" + d5.frontend + "_" + d5.backend + "_" + str(int(d5.bw))
                #################################################################
                toa_flags['pta'] = "InPTA"
                if (epoch_band3.intday()-epoch_band5.intday() > 1):
                    print("epochs are more than 1 day apart")
                    print(ok)
                ar_mjd = epoch_band3.intday()
                ar_freq3 = np.mean(d3.freqs)
                ar_bec3 = d3.beconfig
                ar_freq5 = np.mean(d5.freqs)
                ar_bec5 = d5.beconfig
                if (ar_mjd < 58600.):
                    toa_flags['group'] = 'GM_GWB_'+ar_bec3[1]+'_'+ar_bec3[2]+'_'+ar_bec3[0]+'_pre36'+'_&_'+'GM_GWB_'+ar_bec5[1]+'_'+ar_bec5[2]+'_'+ar_bec5[0]+'_pre36'
                    toa_flags['sys'] = 'GM_GWB_'+ar_bec3[1]+'_'+ar_bec3[2]+'_'+ar_bec3[0]+'_&_'+'GM_GWB_'+ar_bec5[1]+'_'+ar_bec5[2]+'_'+ar_bec5[0]
                if (ar_mjd > 58600.):
                    toa_flags['group'] = 'GM_GWB_'+ar_bec3[1]+'_'+ar_bec3[2]+'_'+ar_bec3[0]+'_post36'+'_&_'+'GM_GWB_'+ar_bec5[1]+'_'+ar_bec5[2]+'_'+ar_bec5[0]+'_post36'
                    toa_flags['sys'] = 'GM_GWB_'+ar_bec3[1]+'_'+ar_bec3[2]+'_'+ar_bec3[0]+'_&_'+'GM_GWB_'+ar_bec5[1]+'_'+ar_bec5[2]+'_'+ar_bec5[0]
                
                if (300. < ar_freq3 < 500.):
                    bnda = '3'
                if (525. < ar_freq3 < 1000.):
                    bnda = '4'
                if (1260. < ar_freq3 < 1460.):
                    bnda = '5'
                if (300. < ar_freq5 < 500.):
                    bndb = '3'
                if (525. < ar_freq5 < 1000.):
                    bndb = '4'
                if (1260. < ar_freq5 < 1460.):
                    bndb = '5'
                toa_flags['bandno'] = bnda + '_&_' + bndb

                if (ar_mjd > 58230.):
                    toa_flags['cycle'] = 'post34'
                    # ext_flag = '-cycle_post34'
                if (ar_mjd < 58230.):
                    toa_flags['cycle'] = 'pre34'
                    # ext_flag = '-cycle_pre34'
                # print(ok)
                #################################################################
                toa_flags['nbin_band3'] = nbin_band3
                toa_flags['nch_band3'] = nchan_band3
                toa_flags['nchx_band3'] = len(freqsx_band3)
                toa_flags['bw_band3'] = freqsx_band3.max() - freqsx_band3.min()
                toa_flags['chbw_band3'] = old_div(abs(d3.bw), nchan_band3)
                toa_flags['subint'] = isub
                toa_flags['tobs_band3'] = d3.subtimes[isub]
                toa_flags['fratio_band3'] = old_div(freqsx_band3.max(), freqsx_band3.min())
                toa_flags['tmplt_band3'] = self.modelfile_band3
                toa_flags['nbin_band5'] = nbin_band5
                toa_flags['nch_band5'] = nchan_band5
                toa_flags['nchx_band5'] = len(freqsx_band5)
                toa_flags['bw_band5'] = freqsx_band5.max() - freqsx_band5.min()
                toa_flags['chbw_band5'] = old_div(abs(d5.bw), nchan_band5)
                toa_flags['tobs_band5'] = d5.subtimes[isub]
                toa_flags['fratio_band5'] = old_div(freqsx_band5.max(), freqsx_band5.min())
                toa_flags['tmplt_band5'] = self.modelfile_band5
                toa_flags['snr'] = results.snr
                if (nu_ref_DM is not None and np.all(fit_flags[:2])):
                    toa_flags['phi_DM_cov'] = \
                        results.covariance_matrix[0, 1]
                toa_flags['gof'] = results.red_chi2
                if print_phase:
                    #########################################################################################################
                    toa_flags['phi'] = results.phi
                    toa_flags['phi_err'] = results.phi_err
                    toa_flags['phi_inf'] = results.phi_inf
                    toa_flags['phi_inf_err'] = results.phi_inf_err
                    toa_flags['DM'] = results.DM
                    toa_flags['DM_err'] = results.DM_err
                    toa_flags['TOA'] = "%d" % (results.TOA.intday()) + ("%.15f" % (results.TOA.fracday()))[1:]
                    toa_flags['TOA_err'] = results.TOA_err
                    toa_flags['TOA_inf'] = "%d" % (results.TOA_inf.intday()) + ("%.15f" % (results.TOA_inf.fracday()))[1:]
                    toa_flags['TOA_inf_err'] = results.TOA_inf_err
                    toa_flags['nu_DM'] = results.nu_DM
                    toa_flags['phi_DM_cov'] = \
                        results.covariance_matrix[0, 1]
                    toa_flags['epoch_b3'] = "%d" % (epoch_band3.intday()) + ("%.15f" % (epoch_band3.fracday()))[1:]
                    toa_flags['epoch_b5'] = "%d" % (epoch_band5.intday()) + ("%.15f" % (epoch_band5.fracday()))[1:]
                    toa_flags['be_delay_b3'] = d3.backend_delay
                    toa_flags['be_delay_b5'] = d5.backend_delay
                    toa_flags['P_b3'] = P_band3
                    toa_flags['P_b5'] = P_band5
                    toa_flags['phase_offset'] = phase_offset
                if bayes:
                        toa_flags['log_z'] = results.log_z
                        for key in results.gt:
                            toa_flags[key] = results.gt[key]
                if bayes and mlan:
                    toa_flags['technique'] = "MLAN_bayes_B35"
                if not bayes and mlan:
                    toa_flags['technique'] = "MLAN_freq_B35"
                if bayes and not mlan:
                    toa_flags['technique'] = "MLA_bayes_B35"
                if not bayes and not mlan:
                    toa_flags['technique'] = "MLA_freq_B35"
                    ###########################################################################################################
                if print_flux:
                    toa_flags['flux'] = fluxes[isub]
                    # consistent with pat / psrflux
                    toa_flags['flux_err'] = flux_errs[isub]
                    #toa_flags['fluxe'] = flux_errs[isub]
                    toa_flags['flux_ref_freq'] = flux_freqs[isub]
                if print_parangle:
                    toa_flags['par_angle'] = d3.parallactic_angles[isub]
                for k, v in addtnl_toa_flags.items():
                    toa_flags[k] = v
                
                self.TOA_list.append(TOA(datafile_band3, datafile_band5, results.nu_DM, results.TOA,
                                         results.TOA_err, d3.telescope, d5.telescope, 
                                         d3.telescope_code, d5.telescope_code, results.DM,
                                         results.DM_err, toa_flags))
                itoa += 1

            DeltaDMs = DMs - DM0
            # The below returns the weighted mean and the sum of the weights,
            # but needs to do better in the case of small-error outliers from
            # RFI, etc.  Also, last TOA may mess things up...use median...?
            
            if np.all(DM_errs[d3.ok_isubs]):
                DM_weights = DM_errs[d3.ok_isubs] ** -2
            else:
                DM_weights = np.ones(len(DM_errs[d3.ok_isubs]))
            DeltaDM_mean, DeltaDM_var = np.average(DeltaDMs[d3.ok_isubs],
                                                   weights=DM_weights, returned=True)
            DeltaDM_var = DeltaDM_var ** -1
            if len(d3.ok_isubs) > 1:
                # The below multiply by the red. chi-squared to inflate the
                # errors.
                DeltaDM_var *= old_div(np.sum(
                    ((DeltaDMs[d3.ok_isubs] - DeltaDM_mean) ** 2) * DM_weights), (len(DeltaDMs[d3.ok_isubs]) - 1))
            DeltaDM_err = DeltaDM_var ** 0.5
            
            self.order_band3.append(datafile_band3)
            self.order_band5.append(datafile_band5)
            self.obs_band3.append(obs_band3)
            self.obs_band5.append(obs_band5)
            self.doppler_fs_band3.append(d3.doppler_factors)
            self.doppler_fs_band5.append(d5.doppler_factors)
            self.nu0s_band3.append(d3.nu0)
            self.nu0s_band5.append(d5.nu0)
            self.nu_fits.append(nu_fits)
            self.nu_refs.append(nu_refs)
            self.ok_isubs.append(d3.ok_isubs)
            self.epochs_band3.append(d3.epochs)
            self.epochs_band5.append(d5.epochs)
            self.MJDs.append(MJDs_band3)
# =============================================================================
#             self.MJDs_band5.append(MJDs_band5)
# =============================================================================
            self.Ps.append(d3.Ps)
# =============================================================================
#             self.Ps_band5.append(d5.Ps)
# =============================================================================
            self.phis.append(phis)  # NB: phis are w.r.t. nu_ref_DM
            self.phi_errs.append(phi_errs)
            self.TOAs.append(TOAs)  # NB: TOAs are w.r.t. nu_ref_DM
            self.TOA_errs.append(TOA_errs)
            self.DM0s.append(DM0)
            self.DMs.append(DMs)
            self.DM_errs.append(DM_errs)
            self.DeltaDM_means.append(DeltaDM_mean)
            self.DeltaDM_errs.append(DeltaDM_err)
            self.GMs.append(GMs)
            self.GM_errs.append(GM_errs)
            self.taus.append(taus)  # NB: taus are w.r.t. nu_ref_tau
            self.tau_errs.append(tau_errs)
            self.alphas.append(alphas)
            self.alpha_errs.append(alpha_errs)
            self.scales.append(scales)
            self.scale_errs.append(scale_errs)
            self.snrs.append(snrs)
            self.channel_snrs.append(channel_snrs)
            self.profile_fluxes.append(profile_fluxes)
            self.profile_flux_errs.append(profile_flux_errs)
            self.fluxes.append(fluxes)
            self.flux_errs.append(flux_errs)
            self.flux_freqs.append(flux_freqs)
            self.covariances.append(covariances)
            self.red_chi2s.append(red_chi2s)
            self.nfevals.append(nfevals)
            self.rcs.append(rcs)
            self.fit_durations.append(fit_duration)
            if not quiet:
                print("--------------------------")
                print(datafile_band3)
                print(datafile_band5)
                print("~%.4f sec/TOA" % (old_div(fit_duration, len(d3.ok_isubs))))
                print("Med. TOA error is %.3f us" % (np.median(
                    phi_errs[d3.ok_isubs]) * d3.Ps.mean() * 1e6))
            if show_plot:
                stop = time.time()
                tot_duration += stop - start
                for isub in d3.ok_isubs:
                    if show_plot=="save":
                        savefig_band3 = datafile_band3 + '.%d3.pptoas.png' % isub
                        savefig_band5 = datafile_band5 + '.%d5.pptoas.png' % isub
                    else:
                        savefig_band3 = False
                        savefig_band5 = False
                    self.show_fit_band3(datafile_band3, isub, savefig=savefig_band3)
                    self.show_fit_band5(datafile_band5, isub, savefig=savefig_band5, phase_offset=phase_offset)
                start = time.time()
            if not show_plot:
                tot_duration = time.time() - start
            if not quiet and len(self.ok_isubs):
                print("--------------------------")
                print("Total time: %.2f sec, ~%.4f sec/TOA" % (tot_duration,
                                                               old_div(tot_duration,
                                                                       (np.array(
                                                                           list(map(len, self.ok_isubs))).sum()))))

    def get_narrowband_TOAs(self, datafile=None, tscrunch=False,
            fit_scat=False, log10_tau=True, scat_guess=None, print_phase=False,
            print_flux=False, print_parangle=False,
            add_instrumental_response=False, addtnl_toa_flags={},
            method='trust-ncg', bounds=None, show_plot=False, quiet=None):
        """
        Measure narrowband TOAs using internal algorithm.

        datafile defaults to self.datafiles, otherwise it is a single
            PSRCHIVE archive name
        tscrunch=True tscrunches archive before fitting (i.e. make one set of
            measurements per archive)
        fit_scat=True will fit the scattering timescale for each TOA. [NOT YET
            IMPLEMENTED]
        log10_tau=True does the scattering fit with log10(scattering timescale)
            as the parameter.
        scat_guess can be a list of three numbers: a guess of the scattering
            timescale tau [s], its reference frequency [MHz], and a guess of
            the scattering index alpha.  Will be used for all archives;
            supercedes other initial values.
        print_phase=True will print the fitted parameter phi and its
            uncertainty on the TOA line with the flags -phs and -phs_err.
        print_flux=True will print an estimate of the flux density and its
            uncertainty on the TOA line.
        print_parangle=True will print the parallactic angle on the TOA line.
        add_instrumental_response=True will account for the instrumental
            response according to the dictionary instrumental_response_dict.
        addtnl_toa_flags are pairs making up TOA flags to be written uniformly
            to all IPTA-formatted TOAs. e.g. ('pta','NANOGrav','version',0.1)
        method is the scipy.optimize.minimize method; currently can be 'TNC',
            'Newton-CG', or 'trust-cng', which are all Newton
            Conjugate-Gradient algorithms.
        bounds is a list of two 2-tuples, giving the lower and upper bounds on
            the phase and tau parameters,
            respectively.  NB: this is only used if method=='TNC'.
        show_plot=True will show a plot of the fitted model, data, and
            residuals at the end of the fitting. [NOT YET IMPLEMENTED]
        quiet=True suppresses output.
        """
        if quiet is None: quiet = self.quiet
        already_warned = False
        warning_message = \
                "You are using an experimental functionality of pptoas!"
        self.nfit = 1
        if fit_scat: self.nfit += 2
        self.fit_phi = True
        self.fit_tau = fit_scat
        self.fit_flags = [int(self.fit_phi), int(self.fit_tau)]
        self.log10_tau = log10_tau
        if not fit_scat:
            self.log10_tau = log10_tau = False
        if True:#fit_scat or self.fit_tau:
            print(warning_message)
            already_warned = True
        self.scat_guess = scat_guess
        start = time.time()
        tot_duration = 0.0
        if datafile is None:
            datafiles = self.datafiles
        else:
            datafiles = [datafile]
        self.tscrunch = tscrunch
        self.add_instrumental_response = add_instrumental_response
        for iarch, datafile in enumerate(datafiles):
            fit_duration = 0.0
            #Load data
            try:
                data = load_data(datafile, dedisperse=False,
                        dededisperse=False, tscrunch=tscrunch,
                        pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                        flux_prof=False, refresh_arch=False, return_arch=False,
                        quiet=quiet)
                if data.dmc:
                    if not quiet:
                        print("%s is dedispersed (dmc = 1).  Reloading it."%\
                                datafile)
                    #continue
                    data = load_data(datafile, dedisperse=False,
                            dededisperse=True, tscrunch=tscrunch,
                            pscrunch=True, fscrunch=False,
                            rm_baseline=rm_baseline, flux_prof=False,
                            refresh_arch=False, return_arch=False, quiet=quiet)
                if not len(data.ok_isubs):
                    if not quiet:
                        print("No subints to fit for %s.  Skipping it."%\
                                datafile)
                    continue
                else: self.ok_idatafiles.append(iarch)
            except RuntimeError:
                if not quiet:
                    print("Cannot load_data(%s).  Skipping it."%datafile)
                continue
            # Unpack the data dictionary into the local namespace; see load_data
            # for dictionary keys.
            # for key in data.keys():
            #     exec(key + " = data['" + key + "']")
            # BWM: Python 3 no longer supports updating local namespace via exec statement/function
            # Also, since we're setting local variables, it's not safe to just update the object __dict__
            # So the best option is to instead load the dictionary into a SimpleNamespace and
            # access key/values as you would an object attribute
            d = SimpleNamespace(**data)
            nsub = d.nsub
            nchan = d.nchan
            nbin = d.nbin
            npol = d.npol
            if d.source is None: d.source = "noname"
            #Observation info
            obs = DataBunch(telescope=d.telescope, backend=d.backend,
                    frontend=d.frontend)
            phis = np.zeros([nsub, nchan], dtype=np.double)
            phi_errs = np.zeros([nsub, nchan], dtype=np.double)
            TOAs = np.zeros([nsub, nchan], dtype="object")
            TOA_errs = np.zeros([nsub, nchan], dtype="object")
            taus = np.zeros([nsub, nchan], dtype=np.float64)
            tau_errs = np.zeros([nsub, nchan], dtype=np.float64)
            scales = np.zeros([nsub, nchan], dtype=np.float64)
            scale_errs = np.zeros([nsub, nchan], dtype=np.float64)
            channel_snrs = np.zeros([nsub, nchan], dtype=np.float64)
            profile_fluxes = np.zeros([nsub, nchan], dtype=np.float64)
            profile_flux_errs = np.zeros([nsub, nchan], dtype=np.float64)
            channel_red_chi2s = np.zeros([nsub, nchan], dtype=np.float64)
            covariances = np.zeros([nsub, nchan, self.nfit, self.nfit],
                    dtype=np.float64)
            nfevals = np.zeros([nsub, nchan], dtype="int")
            rcs = np.zeros([nsub, nchan], dtype="int")
            #PSRCHIVE epochs are *midpoint* of the integration
            MJDs = np.array([d.epochs[isub].in_days() \
                    for isub in range(nsub)], dtype=np.double)
            if self.is_FITS_model:
                if not already_warned:
                    print(warning_message)
                    already_warned = True
                model_data = load_data(self.modelfile, dedisperse=False,
                    dededisperse=False, tscrunch=True, pscrunch=True,
                    fscrunch=False, rm_baseline=True, flux_prof=False,
                    #fscrunch=False, rm_baseline=False, flux_prof=False,
                    refresh_arch=False, return_arch=False, quiet=True)
                model = (model_data.masks * model_data.subints)[0,0]
                if model_data.nbin != nbin:
                    print("Model nbin %d != data nbin %d for %s; skipping it."\
                            %(model_data.nbin, nbin, datafile))
                    continue
                if model_data.nchan == 1:
                    model = np.tile(model[0], len(freqs[0])).reshape(
                            len(freqs[0]), nbin)
                if model_data.nchan != nchan:
                    print("Model nchan %d != data nchan %d for %s; skipping it."%(model_data.nchan, nchan, datafile))
                    continue
            icount = 1
            for isub in d.ok_isubs:
                if self.is_FITS_model and \
                        np.any(model_data.freqs[0] - d.freqs[isub]): # tscrunched
                            print("Frequency mismatch between template and data!")
                sub_id = datafile + "_%d"%isub
                epoch = d.epochs[isub]
                MJD = MJDs[isub]
                P = d.Ps[isub]
                if not self.is_FITS_model:
                    #Read model
                    try:
                        if not fit_scat:
                            self.model_name, self.ngauss, model = read_model(
                                    self.modelfile, d.phases, d.freqs[isub],
                                    d.Ps[isub],
                                    quiet=bool(quiet+(icount-1)))
                        else:
                            self.model_name, self.ngauss, full_model = \
                                    read_model(self.modelfile, d.phases,
                                            d.freqs[isub], d.Ps[isub],
                                            quiet=bool(quiet+(icount-1)))
                            self.model_name, self.model_code, \
                                    self.model_nu_ref, self.ngauss, \
                                    self.gparams, model_fit_flags, self.alpha,\
                                    model_fit_alpha = read_model(
                                            self.modelfile,
                                            quiet=bool(quiet+(icount-1)))
                            unscat_params = np.copy(self.gparams)
                            unscat_params[1] = 0.0
                            model = unscat_model = gen_gaussian_portrait(
                                    self.model_code, unscat_params, 0.0,
                                    d.phases, d.freqs[isub], self.model_nu_ref)
                    except (UnboundLocalError,UnicodeDecodeError):
                        self.model_name, model = read_spline_model(
                                self.modelfile, d.freqs[isub], nbin,
                                quiet=True) #bool(quiet+(icount-1)))
                freqsx = d.freqs[isub,d.ok_ichans[isub]]
                weightsx = d.weights[isub,d.ok_ichans[isub]]
                portx = d.subints[isub,0,d.ok_ichans[isub]]
                modelx = model[d.ok_ichans[isub]]
                if add_instrumental_response and \
                        (self.ird['DM'] or len(self.ird['wids'])):
                            inst_resp_port_FT = instrumental_response_port_FT(
                                    nbin, freqsx, self.ird['DM'], P,
                                    self.ird['wids'], self.ird['irf_types'])
                            modelx = fft.irfft(inst_resp_port_FT * \
                                    fft.rfft(modelx, axis=-1), axis=-1)
                #NB: Time-domain uncertainties below
                errs = d.noise_stds[isub,0,d.ok_ichans[isub]]

                ###################
                # INITIAL GUESSES #
                ###################
                tau_guess = 0.0
                alpha_guess = 0.0
                if fit_scat:
                    if self.scat_guess is not None:
                        tau_guess_s,tau_guess_ref,alpha_guess = self.scat_guess
                        tau_guess = (tau_guess_s / P) * \
                                (nu_fit_tau / tau_guess_ref)**alpha_guess
                    else:
                        if hasattr(self, 'alpha'): alpha_guess = self.alpha
                        else: alpha_guess = scattering_alpha
                        if hasattr(self, 'gparams'):
                            tau_guess = (self.gparams[1] / P) * \
                                    (nu_fit_tau/self.model_nu_ref)**alpha_guess
                        else:
                            tau_guess = 0.0  # nbin**-1?
                            #tau_guess = guess_tau(...)
                    model_prof_scat = fft.irfft(scattering_portrait_FT(
                        np.array([scattering_times(tau_guess, alpha_guess,
                            nu_fit_tau, nu_fit_tau)]), nbin)[0] * fft.rfft(
                                modelx.mean(axis=0)))
                    phi_guess = fit_phase_shift(rot_prof,
                            model_prof_scat, Ns=100).phase
                    if self.log10_tau:
                        if tau_guess == 0.0: tau_guess = nbin**-1
                        tau_guess = np.log10(tau_guess)
                else:
                    #NB: Ns should be larger than nbin for very low S/N data,
                    #especially in the case of noisy models...
                    #phi_guess = fit_phase_shift(rot_prof, modelx.mean(axis=0),
                    #        Ns=100).phase
                    phi_guess = 0.0
                #Need a status bar?
                param_guesses = [phi_guess, tau_guess]
                if bounds is None and method == 'TNC':
                    phi_bounds = (None, None)
                    if not self.log10_tau: tau_bounds = (0.0, None)
                    else: tau_bounds = (np.log10((10*nbin)**-1), None)
                    bounds = [phi_bounds, tau_bounds]

                ###########
                # THE FIT #
                ###########
                fit_flags = list(np.copy(self.fit_flags))
                for ichanx in range(len(d.ok_ichans[isub])):
                    ichan = d.ok_ichans[isub][ichanx]
                    if self.is_FITS_model:
                        if model_data.weights[isub,ichan] == 0: continue
                    prof = portx[ichanx]
                    model_prof = modelx[ichanx]
                    err = errs[ichanx]
                    #results = fit_phase_shift_scat(prof, model_prof,
                    #        param_guesses, P, err, fit_flags, bounds,
                    #        self.log10_tau, option=0, sub_id=sub_id,
                    #        method=method, is_toa=True, quiet=quiet)
                    results = fit_phase_shift(prof, model_prof, err,
                            bounds=[-0.5,0.5], Ns=100)
                    results.tau = results.tau_err = 0.0
                    fit_duration += results.duration

                    ####################
                    #  CALCULATE  TOA  #
                    ####################
                    results.TOA = epoch + pr.MJD(
                            #((results.phi * P) + backend_delay) /
                            ((results.phase * P) + d.backend_delay) /
                            (3600 * 24.))
                    #results.TOA_err = results.phi_err * P * 1e6 # [us]
                    results.TOA_err = results.phase_err * P * 1e6 # [us]

                    #################
                    # ESTIMATE FLUX #
                    #################
                    if print_flux:
                        if results.tau != 0.0:
                            if self.log10_tau: tau = 10**results.tau
                            else: tau = results.tau
                            alpha = results.alpha
                            scat_model = fft.irfft(scattering_profile_FT(tau,
                                data.nbin) * fft.rfft(model_prof))
                        else: scat_model = np.copy(model_prof)
                        scat_model_mean = scat_model.mean()
                        profile_fluxes[isub, ichan] = scat_model_mean * \
                                results.scale
                        profile_flux_errs[isub, ichan] = \
                                abs(scat_model_mean) * results.scale_err

                    phis[isub, ichan] = results.phase
                    phi_errs[isub, ichan] = results.phase_err
                    TOAs[isub, ichan] = results.TOA
                    TOA_errs[isub, ichan] = results.TOA_err
                    taus[isub, ichan] = results.tau
                    tau_errs[isub, ichan] = results.tau_err
                    #nfevals[isub, ichan] = results.nfeval
                    #rcs[isub, ichan] = results.return_code
                    scales[isub, ichan] = results.scale
                    scale_errs[isub, ichan] = results.scale_err
                    channel_snrs[isub, ichan] = results.snr
                    #try:
                    #    covariances[isub, ichan] = results.covariance_matrix
                    #except ValueError:
                    #    for ii,ifit in enumerate(np.where(fit_flags)[0]):
                    #        for jj,jfit in enumerate(np.where(fit_flags)[0]):
                    #            covariances[isub, ichan][ifit,jfit] = \
                    #                    results.covariance_matrix[ii,jj]
                    channel_red_chi2s[isub] = results.red_chi2
                    #Compile useful TOA flags
                    # Add doppler_factor?
                    toa_flags = {}
                    if fit_flags[1]:
                        if self.log10_tau:
                            df = doppler_factors[isub]
                            toa_flags['scat_time'] = 10**results.tau * \
                                    P / df * 1e6  # usec, w/ df
                            toa_flags['log10_scat_time'] = results.tau + \
                                    np.log10(P / df)  # w/ df
                            toa_flags['log10_scat_time_err'] = results.tau_err
                        else:
                            toa_flags['scat_time'] = results.tau * P / df * 1e6
                                                     # usec, w/ df
                            toa_flags['scat_time_err'] = results.tau_err * \
                                    P / df * 1e6  # usec, w/ df
                        toa_flags['phi_tau_cov'] = \
                                results.covariance_matrix[0,1]
                    toa_flags['be'] = d.backend
                    toa_flags['fe'] = d.frontend
                    toa_flags['f'] = d.frontend + "_" + d.backend
                    toa_flags['nbin'] = nbin
                    #toa_flags['nch'] = nchan
                    toa_flags['bw'] = abs(d.bw) / nchan
                    toa_flags['subint'] = isub
                    toa_flags['chan'] = ichan
                    toa_flags['tobs'] = d.subtimes[isub]
                    toa_flags['tmplt'] = self.modelfile
                    toa_flags['snr'] = results.snr
                    toa_flags['gof'] = results.red_chi2
                    if print_phase:
                        toa_flags['phs'] = results.phi
                        toa_flags['phs_err'] = results.phi_err
                    if print_flux:
                        toa_flags['flux'] = fluxes[isub]
                        # consistent with pat / psrflux
                        toa_flags['flux_err'] = flux_errs[isub]
                        #toa_flags['fluxe'] = flux_errs[isub]
                    if print_parangle:
                        toa_flags['par_angle'] = d.parallactic_angles[isub]
                    for k,v in iter(addtnl_toa_flags.items()):
                        toa_flags[k] = v
                    self.TOA_list.append(TOA(datafile, d.freqs[isub,ichan],
                        results.TOA, results.TOA_err, d.telescope,
                        d.telescope_code, None, None, toa_flags))
                icount += 1

            self.order.append(datafile)
            self.obs.append(obs)
            self.doppler_fs.append(d.doppler_factors)
            self.ok_isubs.append(d.ok_isubs)
            self.epochs.append(d.epochs)
            self.MJDs.append(MJDs)
            self.Ps.append(d.Ps)
            self.phis.append(phis)
            self.phi_errs.append(phi_errs)
            self.TOAs.append(TOAs)
            self.TOA_errs.append(TOA_errs)
            self.taus.append(taus)
            self.tau_errs.append(tau_errs)
            self.scales.append(scales)
            self.scale_errs.append(scale_errs)
            self.channel_snrs.append(channel_snrs)
            self.profile_fluxes.append(profile_fluxes)
            self.profile_flux_errs.append(profile_flux_errs)
            self.covariances.append(covariances)
            self.channel_red_chi2s.append(channel_red_chi2s)
            self.nfevals.append(nfevals)
            self.rcs.append(rcs)
            self.fit_durations.append(fit_duration)
            if not quiet:
                print("--------------------------")
                print(datafile)
                print("~%.4f sec/TOA"%(fit_duration / len(self.TOA_list)))
                print("Med. TOA error is %.3f us"%(np.median(
                    phi_errs[d.ok_isubs]) * d.Ps.mean() * 1e6))
            if show_plot:
                pass
                #stop = time.time()
                #tot_duration += stop - start
                #for isub in ok_isubs:
                #    self.show_fit(datafile, isub)
                #start = time.time()
        if not show_plot:
            tot_duration = time.time() - start
        if not quiet and len(self.ok_isubs):
            print("--------------------------")
            print("Total time: %.2f sec, ~%.4f sec/TOA"%(tot_duration,
                    tot_duration / len(self.TOA_list)))

    def get_psrchive_TOAs(self, datafile=None, tscrunch=False, algorithm='PGS',
            toa_format='tempo2', flags='IPTA', attributes=['chan','subint'],
            quiet=False):
        """
        Measure narrowband TOAs using psrchive.

        datafile defaults to self.datafiles, otherwise it is a single
            PSRCHIVE archive name.
        tscrunch=True tscrunches archive before fitting (i.e. make one set of
            measurements per archive).
        algorithm is one of the three-letter codes specifying the shift
            algorithm; see help for 'pat' program, arguments to -A [default =
            'PGS'].
        toa_format and flags are the TOA format and added metadata flags; see
            help for 'pat' program, arguments to -f [default = tempo2 format
            with IPTA flags].
        attributes is a list containing strings that indicate additional TOA
            flags to be included [default = ['chan', 'subint']].
        """
        if quiet is None: quiet = self.quiet
        already_warned = False
        warning_message = \
                "You are using an experimental functionality of pptoas!"
        if True:
            print(warning_message)
            already_warned = True
        self.psrchive_toas = []
        arrtim = pr.ArrivalTime()
        arrtim.set_shift_estimator(algorithm)
        arrtim.set_format(toa_format)
        arrtim.set_format_flags(flags)
        arrtim.set_attributes(attributes)
        if datafile is None:
            datafiles = self.datafiles
        else:
            datafiles = [datafile]
        if self.is_FITS_model:
            model_arch = pr.Archive_load(self.modelfile)
            model_arch.pscrunch()
            arrtim.set_standard(model_arch)
        for iarch, datafile in enumerate(datafiles):
            arch = pr.Archive_load(datafile)
            arch.pscrunch()
            if tscrunch: arch.tscrunch()
            arrtim.set_observation(arch)
            if not self.is_FITS_model:
                nsub,npol,nchan,nbin = arch.get_nsubint(), arch.get_npol(), \
                        arch.get_nchan(),arch.get_nbin()
                freqs = np.array([[sub.get_centre_frequency(ichan) for ichan \
                        in range(nchan)] for sub in arch])
                for isub in range(nsub):
                    #Read model
                    try:
                        Ps = np.array([sub.get_folding_period() for sub in arch],
                                dtype=np.double)
                        phases = get_bin_centers(nbin, lo=0.0, hi=1.0)
                        self.model_name, self.ngauss, model = read_model(
                                self.modelfile, phases, freqs[isub], Ps[isub],
                                quiet=True)
                    except (UnboundLocalError,UnicodeDecodeError):
                        self.model_name, model = read_spline_model(
                                self.modelfile, freqs[isub], nbin, quiet=True)
                    model_arch = arch.clone()
                    model_arch.set_filename(self.modelfile)
                    model_arch.tscrunch()
                    for model_isub in range(1):#nsub):
                        sub = model_arch.get_Integration(model_isub)
                        for ipol in range(npol):
                            for ichan in range(nchan):
                                prof = sub.get_Profile(ipol,ichan)
                                prof.get_amps()[:] = model[ichan]
                                sub.set_weight(ichan, 1.0)
                    arrtim.set_standard(model_arch)
            self.psrchive_toas.append(arrtim.get_toas())

    def get_channels_to_zap(self, SNR_threshold=8.0, rchi2_threshold=1.3,
                            iterate=True, show=False):
        """
        NB: get_TOAs(...) needs to have been called first.

        SNR_threshold is a signal-to-noise ratio value which is used to flag
            channels for zapping (cf. ppzap.py).  Channels that have a S/N
            values below (SNR_threshold**2 / nchx)**0.5, where nchx is the
            number of channels used in the fit, are added to self.zap_channels.
            NB: only operates if SNR_threshold != 0.0 (individual channels may
            have S/N < 0.0).
        rchi2_threshold is a reduced chi-squared value which is used to flag
            channels for zapping (cf. ppzap.py).  Channels that have a reduced
            chi-squared value above rchi2_threshold are added to
            self.zap_channels.
        iterate=True will iterate over the S/N cut by recalculating the
            effective single-channel S/N threshold and continuing cuts until
            no new channels are cut; this helps to ensure all TOAs will have a
            wideband TOA S/N above SNR_threshold.
        show=True will show the before/after portraits for each subint with
            proposed channels to zap.
        """
        for iarch, ok_idatafile in enumerate(self.ok_idatafiles):
            datafile = self.datafiles[ok_idatafile]
            channel_red_chi2s = []
            zap_channels = []
            for isub in self.ok_isubs[iarch]:
                red_chi2s = []
                bad_ichans = []
                port, model, ok_ichans, freqs, noise_stds = self.show_fit(
                    datafile=datafile, isub=isub, rotate=0.0, show=False,
                    return_fit=True, quiet=True)
                channel_snrs = self.channel_snrs[iarch][isub]
                channel_SNR_threshold = (old_div(SNR_threshold ** 2.0, \
                                                 len(ok_ichans))) ** 0.5
                for ichan, ok_ichan in enumerate(ok_ichans):
                    channel_red_chi2 = get_red_chi2(port[ok_ichan],
                                                    model[ok_ichan], errs=noise_stds[ok_ichan],
                                                    dof=len(port[ok_ichan]) - 2)  # Not sure about dof
                    red_chi2s.append(channel_red_chi2)
                    if channel_red_chi2 > rchi2_threshold:
                        bad_ichans.append(ok_ichan)
                    elif np.isnan(channel_red_chi2):
                        bad_ichans.append(ok_ichan)
                    elif SNR_threshold and \
                            channel_snrs[ok_ichan] < channel_SNR_threshold:
                        bad_ichans.append(ok_ichan)
                    else:
                        pass
                channel_red_chi2s.append(red_chi2s)
                zap_channels.append(bad_ichans)
                if iterate and SNR_threshold and len(bad_ichans):
                    old_len = len(bad_ichans)
                    added_new = True
                    while (added_new and (len(ok_ichans) - len(bad_ichans))):
                        # recalculate threshold after removing channels
                        channel_SNR_threshold = (old_div(SNR_threshold ** 2.0, \
                                                         (len(ok_ichans) - len(bad_ichans)))) ** 0.5
                        for ichan, ok_ichan in enumerate(ok_ichans):
                            if ok_ichan in bad_ichans:
                                continue
                            elif channel_snrs[ok_ichan] < \
                                    channel_SNR_threshold:
                                bad_ichans.append(ok_ichan)
                            else:
                                pass
                        added_new = bool(len(bad_ichans) - old_len)
                        old_len = len(bad_ichans)
                if show and len(bad_ichans):
                    show_portrait(port, get_bin_centers(port.shape[1]),
                                  title="%s, subint: %d\nbad chans: %s" % (datafile,
                                                                           isub, bad_ichans), show=False)
                    port[bad_ichans] *= 0.0
                    show_portrait(port, get_bin_centers(port.shape[1]),
                                  title="%s, subint: %d\nbad chans: %s" % (datafile,
                                                                           isub, bad_ichans), show=True)
            self.channel_red_chi2s.append((channel_red_chi2s))
            self.zap_channels.append((zap_channels))

    def show_subint(self, datafile=None, isub=0, rotate=0.0, quiet=None):
        """
        Plot a phase-frequency portrait of a subintegration.

        datafile is a single PSRCHIVE archive name; defaults to the first one
            listed in self.datafiles.
        isub is the index of the subintegration to be displayed.
        rotate is a phase [rot] specifying the amount to rotate the portrait.
        quiet=True suppresses output.

        To be improved.
        (see show_portrait(...))
        """
        if quiet is None: quiet = self.quiet
        if datafile is None:
            datafile = self.datafiles[0]
        ifile = list(np.array(self.datafiles)[self.ok_idatafiles]).index(
            datafile)
        data = load_data(datafile, dedisperse=True,
                         dededisperse=False, tscrunch=self.tscrunch,
                         # pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                         pscrunch=True, fscrunch=False, rm_baseline=True,
                         flux_prof=False, refresh_arch=False, return_arch=False,
                         quiet=quiet)
        title = "%s ; subint %d" % (datafile, isub)
        port = data.masks[isub, 0] * data.subints[isub, 0]
        if rotate: port = rotate_data(port, rotate)
        show_portrait(port=port, phases=data.phases, freqs=data.freqs[isub],
                      title=title, prof=True, fluxprof=True, rvrsd=bool(data.bw < 0))

    def show_fit(self, datafile=None, isub=0, rotate=0.0, show=True,
                 return_fit=False, savefig=False, quiet=None):
        """
        Plot the fit results from a subintegration.

        datafile is a single PSRCHIVE archive name; defaults to the first one
            listed in self.datafiles.
        isub is the index of the subintegration to be displayed.
        rotate is a phase [rot] specifying the amount to rotate the portrait.
        quiet=True suppresses output.

        To be improved.
        (see show_residual_plot(...))
        """
        if quiet is None: quiet = self.quiet
        if datafile is None:
            datafile = self.datafiles[0]
        ifile = list(np.array(self.datafiles)[self.ok_idatafiles]).index(datafile)
        data = load_data(datafile, dedisperse=False,
                         dededisperse=False, tscrunch=self.tscrunch,
                         # pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                         pscrunch=True, fscrunch=False, rm_baseline=True,
                         flux_prof=False, refresh_arch=False, return_arch=False,
                         quiet=quiet)
        if data.dmc:
            if not quiet:
                print("%s is dedispersed (dmc = 1).  Reloading it." % datafile)
            # continue
            data = load_data(datafile, dedisperse=False,
                             dededisperse=True, tscrunch=self.tscrunch,
                             # pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                             pscrunch=True, fscrunch=False, rm_baseline=True,
                             flux_prof=False, refresh_arch=False, return_arch=False,
                             quiet=quiet)
        phi = self.phis[ifile][isub]
        DM = self.DMs[ifile][isub]
        GM = self.GMs[ifile][isub]
        if self.bary:  # get fitted values
            DM /= self.doppler_fs[ifile][isub]
            GM /= self.doppler_fs[ifile][isub] ** 3
        scales = self.scales[ifile][isub]
        freqs = data.freqs[isub]
        nu_ref_DM, nu_ref_GM, nu_ref_tau = self.nu_refs[ifile][isub]
        P = data.Ps[isub]
        phases = data.phases
        if self.is_FITS_model:
            model_data = load_data(self.modelfile, dedisperse=False,
                                   dededisperse=False, tscrunch=True, pscrunch=True,
                                   fscrunch=False, rm_baseline=True, flux_prof=False,
                                   # fscrunch=False, rm_baseline=False, flux_prof=False,
                                   refresh_arch=False, return_arch=False, quiet=True)
            model = (model_data.masks * model_data.subints)[0, 0]
            if model_data.nchan == 1:
                model = np.tile(model[0], len(freqs)).reshape(len(freqs),
                                                              model_data.nbin)
            model_name = self.modelfile
        else:
            try:
                model_name, ngauss, model = read_model(self.modelfile, phases,
                                                       freqs, data.Ps.mean(), quiet=quiet)
                # freqs, data.Ps[isub], quiet=quiet)     #Track down
                if self.taus[ifile][isub] != 0.0:
                    model_name, model_code, model_nu_ref, ngauss, gparams, \
                    model_fit_flags, model_alpha, model_fit_alpha = \
                        read_model(self.modelfile, quiet=quiet)
                    gparams[1] = 0.0
                    model = gen_gaussian_portrait(model_code, gparams, 0.0,
                                                  phases, freqs, model_nu_ref)
            except:
                model_name, model = read_spline_model(self.modelfile,
                        freqs, data.nbin, quiet=True)
        if self.add_instrumental_response and \
                (self.ird['DM'] or len(self.ird['wids'])):
            inst_resp_port_FT = instrumental_response_port_FT(
                data.nbin, freqs, self.ird['DM'], P,
                self.ird['wids'], self.ird['irf_types'])
            model = fft.irfft(inst_resp_port_FT * fft.rfft(model),
                              axis=-1)
        if self.taus[ifile][isub] != 0.0:
            tau = self.taus[ifile][isub]
            if self.log10_tau: tau = 10 ** tau
            alpha = self.alphas[ifile][isub]
            model = fft.irfft(scattering_portrait_FT(
                scattering_times(tau, alpha, freqs, nu_ref_tau), data.nbin) * \
                              fft.rfft(model, axis=1), axis=1)
        port = rotate_portrait_full(data.subints[isub, 0], phi, DM, GM, freqs,
                                    nu_ref_DM, nu_ref_GM, P)
        if rotate:
            model = rotate_data(model, rotate)
            port = rotate_data(port, rotate)
        port *= data.masks[isub, 0]
        model_scaled = np.transpose(scales * np.transpose(model))
        titles = ("%s\nSubintegration %d" % (datafile, isub),
                  "Fitted Model %s" % (model_name), "Residuals")
        if show:
            show_residual_plot(port=port, model=model_scaled, resids=None,
                               phases=phases, freqs=freqs,
                               noise_stds=data.noise_stds[isub, 0], nfit=2, titles=titles,
                               rvrsd=bool(data.bw < 0), savefig=savefig)
        if return_fit:
            return (port, model_scaled, data.ok_ichans[isub], freqs,
                    data.noise_stds[isub, 0])
        
    def show_fit_band3(self, datafile=None, isub=0, rotate=0.0, show=True,
                    return_fit=False, savefig=False, quiet=None):
            """
            Plot the fit results from a subintegration.

            datafile is a single PSRCHIVE archive name; defaults to the first one
                listed in self.datafiles.
            isub is the index of the subintegration to be displayed.
            rotate is a phase [rot] specifying the amount to rotate the portrait.
            quiet=True suppresses output.

            To be improved.
            (see show_residual_plot(...))
            """
            if quiet is None: quiet = self.quiet
            if datafile is None:
                datafile = self.datafiles_band3[0]
            ifile = list(np.array(self.datafiles_band3)[self.ok_idatafiles_band3]).index(datafile)
            data = load_data(datafile, dedisperse=False,
                            dededisperse=False, tscrunch=self.tscrunch,
                            # pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                            pscrunch=True, fscrunch=False, rm_baseline=True,
                            flux_prof=False, refresh_arch=False, return_arch=False,
                            quiet=quiet)
            if data.dmc:
                if not quiet:
                    print("%s is dedispersed (dmc = 1).  Reloading it." % datafile)
                # continue
                data = load_data(datafile, dedisperse=False,
                                dededisperse=True, tscrunch=self.tscrunch,
                                # pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                                pscrunch=True, fscrunch=False, rm_baseline=True,
                                flux_prof=False, refresh_arch=False, return_arch=False,
                                quiet=quiet)
            phi = self.phis[ifile][isub]
            DM = self.DMs[ifile][isub]
            GM = self.GMs[ifile][isub]
            if self.bary:  # get fitted values
                DM /= self.doppler_fs_band3[ifile][isub]
                GM /= self.doppler_fs_band3[ifile][isub] ** 3
            #scales = self.scales[ifile][isub]
            freqs = data.freqs[isub]
            #This only works if scales are added band3 + band5
            scales = self.scales[ifile][isub][:len(freqs)]
            #
            nu_ref_DM, nu_ref_GM, nu_ref_tau = self.nu_refs[ifile][isub]
            P = data.Ps[isub]
            phases = data.phases
            if self.is_FITS_model_band3:
                model_data = load_data(self.modelfile_band3, dedisperse=False,
                                    dededisperse=False, tscrunch=True, pscrunch=True,
                                    fscrunch=False, rm_baseline=True, flux_prof=False,
                                    # fscrunch=False, rm_baseline=False, flux_prof=False,
                                    refresh_arch=False, return_arch=False, quiet=True)
                model = (model_data.masks * model_data.subints)[0, 0]
                if model_data.nchan == 1:
                    model = np.tile(model[0], len(freqs)).reshape(len(freqs),
                                                                model_data.nbin)
                model_name = self.modelfile_band3
            else:
                try:
                    model_name, ngauss, model = read_model(self.modelfile_band3, phases,
                                                        freqs, data.Ps.mean(), quiet=quiet)
                    # freqs, data.Ps[isub], quiet=quiet)     #Track down
                    if self.taus[ifile][isub] != 0.0:
                        model_name, model_code, model_nu_ref, ngauss, gparams, \
                        model_fit_flags, model_alpha, model_fit_alpha = \
                            read_model(self.modelfile_band3, quiet=quiet)
                        gparams[1] = 0.0
                        model = gen_gaussian_portrait(model_code, gparams, 0.0,
                                                    phases, freqs, model_nu_ref)
                except:
                    model_name, model = read_spline_model(self.modelfile_band3,
                            freqs, data.nbin, quiet=True)
            if self.add_instrumental_response and \
                    (self.ird_band3['DM'] or len(self.ird_band3['wids'])):
                inst_resp_port_FT = instrumental_response_port_FT(
                    data.nbin, freqs, self.ird_band3['DM'], P,
                    self.ird_band3['wids'], self.ird_band3['irf_types'])
                model = fft.irfft(inst_resp_port_FT * fft.rfft(model),
                                axis=-1)
            if self.taus[ifile][isub] != 0.0:
                tau = self.taus[ifile][isub]
                if self.log10_tau: tau = 10 ** tau
                alpha = self.alphas[ifile][isub]
                model = fft.irfft(scattering_portrait_FT(
                    scattering_times(tau, alpha, freqs, nu_ref_tau), data.nbin) * \
                                fft.rfft(model, axis=1), axis=1)
            port = rotate_portrait_full(data.subints[isub, 0], phi, DM, GM, freqs,
                                        nu_ref_DM, nu_ref_GM, P)
            if rotate:
                model = rotate_data(model, rotate)
                port = rotate_data(port, rotate)
            port *= data.masks[isub, 0]
            model_scaled = np.transpose(scales * np.transpose(model))
            titles = ("%s\nSubintegration %d" % (datafile, isub),
                    "Fitted Model %s" % (model_name), "Residuals")
            if show:
                show_residual_plot(port=port, model=model_scaled, resids=None,
                                phases=phases, freqs=freqs,
                                noise_stds=data.noise_stds[isub, 0], nfit=2, titles=titles,
                                rvrsd=bool(data.bw < 0), savefig=savefig)
            if return_fit:
                return (port, model_scaled, data.ok_ichans[isub], freqs,
                        data.noise_stds[isub, 0])

    def show_fit_band5(self, datafile=None, isub=0, rotate=0.0, show=True,
                    return_fit=False, savefig=False, quiet=None, phase_offset=0.0):
            """
            Plot the fit results from a subintegration.

            datafile is a single PSRCHIVE archive name; defaults to the first one
                listed in self.datafiles.
            isub is the index of the subintegration to be displayed.
            rotate is a phase [rot] specifying the amount to rotate the portrait.
            quiet=True suppresses output.

            To be improved.
            (see show_residual_plot(...))
            """
            if quiet is None: quiet = self.quiet
            if datafile is None:
                datafile = self.datafiles_band5[0]
            ifile = list(np.array(self.datafiles_band5)[self.ok_idatafiles_band5]).index(datafile)
            data = load_data(datafile, dedisperse=False,
                            dededisperse=False, tscrunch=self.tscrunch,
                            # pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                            pscrunch=True, fscrunch=False, rm_baseline=True,
                            flux_prof=False, refresh_arch=False, return_arch=False,
                            quiet=quiet)
            if data.dmc:
                if not quiet:
                    print("%s is dedispersed (dmc = 1).  Reloading it." % datafile)
                # continue
                data = load_data(datafile, dedisperse=False,
                                dededisperse=True, tscrunch=self.tscrunch,
                                # pscrunch=True, fscrunch=False, rm_baseline=rm_baseline,
                                pscrunch=True, fscrunch=False, rm_baseline=True,
                                flux_prof=False, refresh_arch=False, return_arch=False,
                                quiet=quiet)
            phi = self.phis[ifile][isub] + phase_offset
            DM = self.DMs[ifile][isub]
            GM = self.GMs[ifile][isub]
            if self.bary:  # get fitted values
                DM /= self.doppler_fs_band5[ifile][isub]
                GM /= self.doppler_fs_band5[ifile][isub] ** 3
            #scales = self.scales[ifile][isub]
            freqs = data.freqs[isub]
            #This only works if scales are added band3 + band5
            scales = self.scales[ifile][isub][-1*len(freqs):]
            #
            nu_ref_DM, nu_ref_GM, nu_ref_tau = self.nu_refs[ifile][isub]
            P = data.Ps[isub]
            phases = data.phases
            if self.is_FITS_model_band5:
                model_data = load_data(self.modelfile_band5, dedisperse=False,
                                    dededisperse=False, tscrunch=True, pscrunch=True,
                                    fscrunch=False, rm_baseline=True, flux_prof=False,
                                    # fscrunch=False, rm_baseline=False, flux_prof=False,
                                    refresh_arch=False, return_arch=False, quiet=True)
                model = (model_data.masks * model_data.subints)[0, 0]
                if model_data.nchan == 1:
                    model = np.tile(model[0], len(freqs)).reshape(len(freqs),
                                                                model_data.nbin)
                model_name = self.modelfile_band5
            else:
                try:
                    model_name, ngauss, model = read_model(self.modelfile_band5, phases,
                                                        freqs, data.Ps.mean(), quiet=quiet)
                    # freqs, data.Ps[isub], quiet=quiet)     #Track down
                    if self.taus[ifile][isub] != 0.0:
                        model_name, model_code, model_nu_ref, ngauss, gparams, \
                        model_fit_flags, model_alpha, model_fit_alpha = \
                            read_model(self.modelfile_band5, quiet=quiet)
                        gparams[1] = 0.0
                        model = gen_gaussian_portrait(model_code, gparams, 0.0,
                                                    phases, freqs, model_nu_ref)
                except:
                    model_name, model = read_spline_model(self.modelfile_band5,
                            freqs, data.nbin, quiet=True)
            if self.add_instrumental_response and \
                    (self.ird_band5['DM'] or len(self.ird_band5['wids'])):
                inst_resp_port_FT = instrumental_response_port_FT(
                    data.nbin, freqs, self.ird_band5['DM'], P,
                    self.ird_band5['wids'], self.ird_band5['irf_types'])
                model = fft.irfft(inst_resp_port_FT * fft.rfft(model),
                                axis=-1)
            if self.taus[ifile][isub] != 0.0:
                tau = self.taus[ifile][isub]
                if self.log10_tau: tau = 10 ** tau
                alpha = self.alphas[ifile][isub]
                model = fft.irfft(scattering_portrait_FT(
                    scattering_times(tau, alpha, freqs, nu_ref_tau), data.nbin) * \
                                fft.rfft(model, axis=1), axis=1)
            port = rotate_portrait_full(data.subints[isub, 0], phi, DM, GM, freqs,
                                        nu_ref_DM, nu_ref_GM, P)
            if rotate:
                model = rotate_data(model, rotate)
                port = rotate_data(port, rotate)
            port *= data.masks[isub, 0]
            model_scaled = np.transpose(scales * np.transpose(model))
            titles = ("%s\nSubintegration %d" % (datafile, isub),
                    "Fitted Model %s" % (model_name), "Residuals")
            if show:
                show_residual_plot(port=port, model=model_scaled, resids=None,
                                phases=phases, freqs=freqs,
                                noise_stds=data.noise_stds[isub, 0], nfit=2, titles=titles,
                                rvrsd=bool(data.bw < 0), savefig=savefig)
            if return_fit:
                return (port, model_scaled, data.ok_ichans[isub], freqs,
                        data.noise_stds[isub, 0])

if __name__ == "__main__":

    from optparse import OptionParser

#Avinash
    usage = "Usage: %prog -a <datafile or metafile of band3> -b <datafile or metafile of band5> -c <modelfile b3> -d <modelfile b5> [options]"
    
    parser = OptionParser(usage)
    # parser.add_option("-h", "--help",
    #                  action="store_true", dest="help", default=False,
    #                  help="Show this help message and exit.")
#Avinash
    parser.add_option("-a", "--datafiles_band3",
                      action="store", metavar="archive", dest="datafiles_band3",
                      default="meta_files/J1909.b3.meta",
                      help="PSRCHIVE archive from which to measure TOAs/DMs, or a metafile listing archive filenames.  \
                              ***Recommended: files should not be dedispersed!*** \
                              (i.e. vap -c dmc <datafile> should return 0)")
    parser.add_option("-b", "--datafiles_band5",
                      action="store", metavar="archive", dest="datafiles_band5",
                      default="meta_files/J1909.b5.meta",
                      help="PSRCHIVE archive from which to measure TOAs/DMs, or a metafile listing archive filenames.  \
                              ***Recommended: files should not be dedispersed!*** \
                              (i.e. vap -c dmc <datafile> should return 0)")
    
    #Avinash                     
    parser.add_option("-c", "--modelfile_band3",
                      action="store", metavar="model", dest="modelfile_band3",
                      default="port_files/J1909.b3.spl",
                      help="Model file from ppgauss.py, ppspline.py, or PSRCHIVE FITS file that either has same channel frequencies, nchan, & nbin as datafile(s), or is a single profile (nchan = 1, with the same nbin) to be interpreted as a constant template.")
    parser.add_option("-d", "--modelfile_band5",
                      action="store", metavar="model", dest="modelfile_band5",
                      default="port_files/J1909.b5.spl",
                      help="Model file from ppgauss.py, ppspline.py, or PSRCHIVE FITS file that either has same channel frequencies, nchan, & nbin as datafile(s), or is a single profile (nchan = 1, with the same nbin) to be interpreted as a constant template.")
    parser.add_option("-o", "--outfile",
                      action="store", metavar="timfile", dest="outfile",
                      default="J1909.pp35_2.tim",
                      help="Name of output .tim file. Will append. [default=stdout]")
    parser.add_option("--narrowband",
                      action="store_true",  dest="narrowband", default=False,
                      help="Make narrowband TOAs instead.")
    parser.add_option("--psrchive",
                      action="store_true",  dest="psrchive", default=False,
                      help="Make narrowband TOAs with PSRCHIVE.")
    parser.add_option("--errfile",
                      action="store", metavar="errfile", dest="errfile",
                      default=None,
                      help="If specified, will write the fitted DM errors to errfile (desirable if using 'Princeton'-like formatted TOAs). Will append.")
    parser.add_option("-T", "--tscrunch",
                      action="store_true", dest="tscrunch", default=False,
                      help="tscrunch archives before measurement (i.e., return only one set of measurements per archive.")
    parser.add_option("-f", "--format",
                      action="store", metavar="format", dest="format",
                      help="Format of output .tim file; either 'princeton' or 'ipta'.  Default is IPTA-like format.")
    parser.add_option("--nu_ref",
                      action="store", metavar="nu_ref", dest="nu_ref_DM",
                      default=None,
                      help="Topocentric frequency [MHz] to which the output TOAs are referenced, i.e. the frequency that has zero delay from a non-zero DM. 'inf' is used as an argument here for infinite frequency, but the default internal behavior follows TEMPO/2 convention and will write 0.0 for infinite-frequency TOAs. [defaults to zero-covariance frequency, recommended]")
    parser.add_option("--DM",
                      action="store", metavar="DM", dest="DM0", default=None,
                      help="Nominal DM [cm**-3 pc] from which to reference offset DM measurements.  If unspecified, will use the DM stored in each archive.")
    parser.add_option("--no_bary",
                      action="store_false", dest="bary", default=True,
                      help="Do not Doppler-correct DMs, GMs, taus, or nu_tau.  Output values are 'topocentric'.")
    parser.add_option("--one_DM",
                      action="store_true", dest="one_DM", default=False,
                      help="Returns single DM value in output .tim file for all subints in the epoch instead of a fitted DM per subint.")
    parser.add_option("--fix_DM",
                      action="store_false", dest="fit_DM", default=True,
                      help="Do not fit for DM. NB: the parfile DM will still be 'barycentered' in the TOA lines unless --no_bary is used!")
    parser.add_option("--fit_dt4",
                      action="store_true", dest="fit_GM", default=False,
                      help="Fit for delays that scale as nu**-4 and return 'GM' parameters s.t. dt4 = Dconst**2 * GM * nu**-4.  GM has units [cm**-6 pc**2 s**-1] and can be related to a discrete cloud causing refractive, geometric delays.")
    parser.add_option("--fit_scat",
                      action="store_true", dest="fit_scat", default=False,
                      help="Fit for scattering timescale and index per TOA.  Can be used with --fix_alpha.")
    parser.add_option("--no_logscat",
                      action="store_false", dest="log10_tau", default=True,
                      help="If using fit_scat, this flag specifies not to fit the log10 of the scattering timescale, but simply the scattering timescale.")
    parser.add_option("--scat_guess",
                      action="store", dest="scat_guess",
                      metavar="tau,freq,alpha",
                      default=None,
                      help="If using fit_scat, manually specify a comma-separated triplet containing an initial guess for the scattering timescale parameter [s], its reference frequency [MHz], and an initial guess for the scattering index.  Will be used for all archives; supercedes other initial values.")
    parser.add_option("--fix_alpha",
                      action="store_true", dest="fix_alpha", default=False,
                      help="Fix the scattering index value to the value specified as scattering_alpha in pplib.py or alpha in the provided .gmodel file.  Only used in combination with --fit_scat.")
    parser.add_option("--nu_tau",
                      action="store", metavar="nu_ref_tau", dest="nu_ref_tau",
                      default=None,
                      help="Frequency [MHz] to which the output scattering times are referenced, i.e. tau(nu) = tau * (nu/nu_ref_tau)**alpha.  If no_bary is True, this frequency is topocentric, otherwise barycentric. [default=nu_zero (zero-covariance frequency, recommended)]")
    parser.add_option("--print_phase",
                      action="store_true", dest="print_phase", default=True,
                      help="Print the fitted phase shift and its uncertainty on the TOA line with the flag -phs")
    parser.add_option("--print_flux",
                      action="store_true", dest="print_flux", default=False,
                      help="Print an estimate of the overall mean flux density and its uncertainty on the TOA line.")
    parser.add_option("--print_parangle",
                      action="store_true", dest="print_parangle",
                      default=False,
                      help="Print the parallactic angle of each subintegration on the TOA line.")
    parser.add_option("--flags",
                      action="store", metavar="flags", dest="toa_flags",
                      default="",
                      help="Pairs making up TOA flags to be written uniformly to all IPTA-formatted TOAs.  e.g., --flags=pta,NANOGrav,version,0.1,key3,val3... etc.")
    parser.add_option("--snr_cut",
                      metavar="S/N", action="store", dest="snr_cutoff",
                      default=0.0,
                      help="Set a S/N cutoff for TOAs written.")
    parser.add_option("--showplot",
                      action="store_true", dest="show_plot", default="save",
                      help="Show a plot of fitted data/model/residuals for each subint.  Good for diagnostic purposes only.")
    parser.add_option("--saveplot",
                      action="store_true", dest="save_plot", default=False,
                      help="Save plots of fitted data/model/residuals for each subint.")
    parser.add_option("--outdir",
                      action="store", dest="outdir", default='./',
                      help="Output directory to store the Bayesian results.")
    parser.add_option("--bayes",
                      action="store_true", dest="bayes", default=False,
                      help="If True, will use bayesian method to estimate ToA and DM")
    parser.add_option("--mlan",
                      action="store_true", dest="mlan", default=False,
                      help="If True, will use MLAN technique to estimate ToA and DM")
    parser.add_option("--pool",
                      action="store", metavar="pool", dest="pool",
                      default=10,
                      help="pool input for nautilus")
    parser.add_option("--quiet",
                      action="store_true", dest="quiet", default=False,
                      help="Only TOAs printed to standard output, if outfile is None.")

    (options, args) = parser.parse_args()
# =============================================================================
#     print("avinash: options",options,type(options))
#     print("avinash: args",args,type(args))
#     print(ok)
# =============================================================================

    if (options.datafiles_band3 is None or options.modelfile_band3 is None or options.datafiles_band5 is None or options.modelfile_band5 is None):
        print("\npptoas.py - simultaneously measure TOAs, DMs, and scattering in broadband data\n")
        parser.print_help()
        print("")
        parser.exit()

    #Avinash
    datafiles_band3 = options.datafiles_band3
    datafiles_band5 = options.datafiles_band5
    modelfile_band3 = options.modelfile_band3
    modelfile_band5 = options.modelfile_band5
    
    outfile = options.outfile
    narrowband = options.narrowband
    psrchive = options.psrchive
    errfile = options.errfile
    tscrunch = options.tscrunch
    format = options.format
    nu_ref_DM = options.nu_ref_DM
    if nu_ref_DM:
        if nu_ref_DM == "inf":
            nu_ref_DM = np.inf
        else:
            nu_ref_DM = np.float64(nu_ref_DM)
        nu_refs = (nu_ref_DM, None)
    else:
        nu_refs = None
    DM0 = options.DM0
    if DM0: DM0 = np.float64(DM0)
    bary = options.bary
    one_DM = options.one_DM
    fit_DM = options.fit_DM
    fit_GM = options.fit_GM
    fit_scat = options.fit_scat
    log10_tau = options.log10_tau
    scat_guess = options.scat_guess
    if scat_guess:
        scat_guess = [s.upper() for s in scat_guess.split(',')]
        scat_guess = list(map(float, scat_guess))
    fix_alpha = options.fix_alpha
    nu_ref_tau = options.nu_ref_tau
    if nu_ref_tau:
        nu_ref_tau = np.float64(nu_ref_tau)
        if nu_ref_DM:
            nu_refs = (nu_ref_DM, nu_ref_tau)
        else:
            nu_refs = (None, nu_ref_tau)
    print_phase = options.print_phase
    print_flux = options.print_flux
    print_parangle = options.print_parangle
    k, v = options.toa_flags.split(',')[::2], options.toa_flags.split(',')[1::2]
    addtnl_toa_flags = dict(list(zip(k, v)))
    snr_cutoff = float(options.snr_cutoff)
    show_plot = options.show_plot
    if options.save_plot: show_plot = "save"
    bayes = options.bayes
    mlan = options.mlan
    pool = options.pool
    quiet = options.quiet
    outdir = options.outdir
    
    if not mlan:
        from mpptoaslib_b35 import *
    if mlan:
        from mpptoaslib_b35_MLAN import *
# =============================================================================
#     print("avinash: datafiles",datafiles_band3)
#     print("avinash: datafiles",datafiles_band5)
#     print("avinash: modelfile",modelfile_band3)
#     print("avinash: modelfile",modelfile_band5)
#     print("avinash: outfile",outfile)
#     print("avinash: narrowband",narrowband)
#     print("avinash: psrchive",psrchive)
#     print("avinash: errfile",errfile)
#     print("avinash: tscrunch",tscrunch)
#     print("avinash: format",format)
#     print("avinash: nu_ref_DM",nu_ref_DM)
#     print("avinash: DM0",DM0)
#     print("avinash: bary",bary)
#     print("avinash: one_DM",one_DM)
#     print("avinash: fit_DM",fit_DM)
#     print("avinash: fit_GM",fit_GM)
#     print("avinash: fit_scat",fit_scat)
#     print("avinash: log10_tau",log10_tau)
#     print("avinash: scat_guess",scat_guess)
#     print("avinash: fix_alpha",fix_alpha)
#     print("avinash: nu_ref_tau",nu_ref_tau)
#     print("avinash: print_phase",print_phase)
#     print("avinash: print_flux",print_flux)
#     print("avinash: print_parangle",print_parangle)
#     print("avinash: k, v",k, v)
#     print("avinash: addtnl_toa_flags",addtnl_toa_flags)
#     print("avinash: snr_cutoff",snr_cutoff)
#     print("avinash: show_plot",show_plot)
#     print("avinash: quiet",quiet)
#     print("avinash: nu_refs",nu_refs)
#     print(ok)
# =============================================================================

    gt = GetTOAs(datafiles_band3=datafiles_band3, datafiles_band5=datafiles_band5, modelfile_band3=modelfile_band3, modelfile_band5=modelfile_band5, quiet=quiet)
    if not narrowband and not psrchive:  # get wideband TOAs
        gt.get_TOAs(datafile_band3=None, datafile_band5=None, tscrunch=tscrunch, nu_refs=nu_refs, DM0=DM0,
                bary=bary, fit_DM=fit_DM, fit_GM=fit_GM, fit_scat=fit_scat,
                log10_tau=log10_tau, scat_guess=scat_guess,
                fix_alpha=fix_alpha, print_phase=print_phase,
                print_flux=print_flux, print_parangle=print_parangle,
                addtnl_toa_flags=addtnl_toa_flags, method='trust-ncg',
                bounds=None, nu_fits=None, show_plot=show_plot, bayes=bayes, quiet=quiet, outdir=outdir)
# =============================================================================
#     elif not psrchive:  # get narrowband TOAs using in-house code
#         gt.get_narrowband_TOAs(datafile=None, tscrunch=tscrunch,
#                 fit_scat=fit_scat, log10_tau=log10_tau, scat_guess=scat_guess,
#                 print_phase=print_phase, print_flux=print_flux,
#                 print_parangle=print_parangle,
#                 addtnl_toa_flags=addtnl_toa_flags, method='trust-ncg',
#                 bounds=None, show_plot=False, quiet=quiet)
#     else:  # get narrowband TOAs with psrchive
#         gt.get_psrchive_TOAs(datafile=None, tscrunch=False, algorithm='PGS',
#                 toa_format='Tempo2', flags='IPTA',
#                 attributes=['chan','subint'])#, print_toas=True, outfile=None)
# =============================================================================
    if not psrchive:
        if format == "princeton":
            gt.write_princeton_TOAs(outfile=outfile, one_DM=one_DM,
                dmerrfile=errfile)
        else:
            if one_DM:
                gt.TOA_one_DM_list = [toa for toa in gt.TOA_list]
                for toa in gt.TOA_one_DM_list:
                    ifile = list(
                            np.array(gt.datafiles)[gt.ok_idatafiles]).index( \
                                    toa_archive)
                    DDM = gt.DeltaDM_means[ifile]
                    DDM_err = gt.DeltaDM_errs[ifile]
                    toa.DM = DDM + gt.DM0s[ifile]
                    toa.DM_error = DDM_err
                    toa.flags['DM_mean'] = True
                write_TOAs(gt.TOA_one_DM_list, inf_is_zero=True,
                        SNR_cutoff=snr_cutoff, outfile=outfile, append=True)
            else:
                print(outfile)
                print(np.shape(gt.TOA_list))
                write_TOAs(gt.TOA_list, inf_is_zero=True, SNR_cutoff=snr_cutoff,
                        outfile=outfile, append=True)
    else:
        if outfile is not None:
            of = open(outfile, 'a')
        for iarch,arch in enumerate(gt.datafiles):
            for itoa,toa in enumerate(gt.psrchive_toas[iarch]):
                if outfile is not None:
                    of.write(gt.psrchive_toas[iarch][itoa]+"\n")
                else:
                    print(gt.psrchive_toas[iarch][itoa])
        if outfile is not None:
            of.close()

