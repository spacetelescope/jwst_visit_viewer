import os
import platform
import subprocess
import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from urllib.parse import quote
import astropy, astropy.wcs, astropy.visualization
import astropy.io.fits as fits
import astropy.coordinates as coords
import astropy.units as u
import pysiaf
import copy


# Visit plotting tools

# A note about interpretation of GSPA. This is subtle, and the implementation varies depending on PPS version.
# OSS passes through the GSPA parameter to spacecraft ACS for slews. ACS does not utilize the V frame in its pointing control.
# - The GSPA parameter is interpreted as the PA of the FGS1 Ideal coordinate system ("FGS ics"), and transformed to the
#   spacecraft J frame via the FGS_to_STA k-constant matrix. In so doing, the GSPA is interpreted as the PA of th FGS Yics
#   angle *at the position of the guide star*
# - Note that the FGS1 Yics differs from the V3 axis by ~ -1.25 degrees
#
# In PPS versions 14.14.1 and below, this was neglected and OPGS provides the computed V3PA at the guide star as the GSPA parameter.
#   (incorrectly / inconsistently, resulting in unexpected attitudes)
# In later versions, PPS applies the compensation and provides the FGS Y ideal PA, which spacecraft ACS can then transform to the J frame
#


PROGRAMS_WITH_SEGMENT_GUIDING = [1410, 1141, 1143, 1148, 1150, 1151, 1153, 1158,  # flight programs, not yet a complete list
                                  710,  741,  743,  # rehearsal programs
                                  646]  # WF Guiding rehearsal program


def get_image_cache_dir():
    cache_dir = os.path.join(os.path.dirname(__file__), "image_cache")
    if not os.path.isdir(cache_dir):
        os.mkdir(cache_dir)
    return cache_dir

def load_all_siafs():
    return {
    'NRC': pysiaf.Siaf("NIRCam"),
    'FGS': pysiaf.Siaf('FGS'),
    'NIS': pysiaf.Siaf('NIRISS'),
    'NRS': pysiaf.Siaf("NIRSpec"),
    'MIRI': pysiaf.Siaf('MIRI')
    }

# let's just load these once and save in a global to reuse as needed, for efficiency
SIAFS = load_all_siafs()

##--- Functions for plotting the visit field of VIEW

def retrieve_2mass_image(visit, ra=None, dec=None, verbose=True, redownload=False, filter='K', fov=0.35):
    """Obtain from Aladin a 2MASS image for the pointing location of a JWST visit

    Uses HIPS2FITS service; see http://alasky.u-strasbg.fr/hips-image-services/hips2fits

    FITS files for the retrieved images are cached for re-use, in a subdirectory
    `image_cache` next to where this code is.

    Parameters
    ----------
    visit : VisitFileContents object
        Representation of some JWST visit file
    filter : string
        Which bandpass filter in 2MASS to get?
    verbose : bool
        more output text
    redownload : bool
        Even if image is already downloaded and cached, ignore that and download from Vizier again.

    """

    if ra==None and dec==None:
        ra = visit.slew.GSRA
        dec = visit.slew.GSDEC

    hips_catalog = f'CDS/P/2MASS/{filter}'  # also try 2MASS/color
    width = 1024
    height = 1024

    visitname = os.path.splitext(os.path.basename(visit.filename))[0]

    if fov!= 0.35:
        # if a non-default FOV is used, save that specially
        img_fn = os.path.join(get_image_cache_dir(),  f'img_2mass_{filter}_{visitname}_fov{fov}.fits')
    else:
        img_fn = os.path.join(get_image_cache_dir(),  f'img_2mass_{filter}_{visitname}.fits')

    if not os.path.exists(img_fn) or redownload:

        # optional / TBD - add PA into this query?
        # rotation_angle=90.0
        url = f'http://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips={quote(hips_catalog)}&width={width}&height={height}&fov={fov}&projection=TAN&coordsys=icrs&ra={ra}&dec={dec}'
        if verbose:
            print(f"Retrieving 2MASS image from Aladin near ra={ra} & dec={dec}...")

        with fits.open(url) as hdu:
            hdu.writeto(img_fn, overwrite=True)
            if verbose:
                print(f"Retrieved 2MASS image from Aladin and saved to {img_fn}")
    hdu = fits.open(img_fn)
    return hdu


def plot_visit_fov(visit, verbose=False, subplotspec=None, use_dss=False, center_visit=None):
    """Make a nice annotated plot of a visit FOV"""


    # let's center the plot on the master chief ray (MCF; between NIRCams) at the science attitude.
    # This works better than centering on the guide star.
    # Use SIAF transforms to infer the MCF pointing in RA, Dec.
    fgs_aperture_name = visit.get_guider_aperture(return_name=True)
    fgs_aperture = SIAFS['FGS'].apertures[fgs_aperture_name]
    attmat = visit.get_attitude_matrix(step='sci')
    fgs_aperture.set_attitude_matrix(attmat)
    mcf_ra, mcf_dec = fgs_aperture.tel_to_sky(0,-468)   # RA, Dec of master chief ray location (between NIRCam A+B)

    if center_visit is not None:
        # Override default center coordinates using some other visit
        print(f"Overriding plot center coordinates; will center based on {center_visit}")
        from . import visitparser
        center_visit_parsed = visitparser.VisitFileContents(center_visit)
        center_fgs_aperture_name = center_visit_parsed.get_guider_aperture(return_name=True)
        center_fgs_aperture = copy.deepcopy(SIAFS['FGS'].apertures[center_fgs_aperture_name])
        center_attmat = center_visit_parsed.get_attitude_matrix(step='sci')
        center_fgs_aperture.set_attitude_matrix(center_attmat)
        mcf_ra, mcf_dec = center_fgs_aperture.tel_to_sky(0, -468)  # RA, Dec of master chief ray location (between NIRCam A+B)

    # Compute RA, Dec, PA of the V1 axis, for comparison to values in SciOps OP delivery report
    v1_ra, v1_dec = fgs_aperture.tel_to_sky(0,0)        # RA, Dec of V1 axis reference location
    v1pv3_ra, v1pv3_dec = fgs_aperture.tel_to_sky(0,1)  # Will use to compute V3PA at V1 axis reference location
    v1c = coords.SkyCoord(v1_ra, v1_dec, unit='deg', frame='icrs')
    v1pv3c = coords.SkyCoord(v1pv3_ra, v1pv3_dec, unit='deg', frame='icrs')
    v3pa_at_v1 = v1c.position_angle(v1pv3c).deg

    # Compute RA, Dec of the J frame axis
    jframe_aperture = SIAFS['FGS'].apertures['J-FRAME']
    jframe_aperture.set_attitude_matrix(attmat)
    j_ra, j_dec = jframe_aperture.idl_to_sky(0,0)        # RA, Dec of V1 axis reference location

    image_visit = center_visit_parsed if center_visit else visit
    if use_dss:
        img_hdu = retrieve_dss_image(image_visit, ra = mcf_ra, dec=mcf_dec, verbose=verbose)
    else:
        img_hdu = retrieve_2mass_image(image_visit, ra = mcf_ra, dec=mcf_dec, verbose=verbose)

    wcs = astropy.wcs.WCS(img_hdu[0].header)

    #-- Setup plot and axes
    if subplotspec is None:
        plt.figure(figsize=(16,9))
        ax = plt.subplot(projection=wcs)
    else:
        ax = plt.subplot(subplotspec,projection=wcs)

    #-- Display the 2MASS image
    norm = astropy.visualization.ImageNormalize(img_hdu[0].data,
                          interval=astropy.visualization.PercentileInterval(99.99),
                      stretch=astropy.visualization.AsinhStretch(a=0.001))
    plt.imshow(img_hdu[0].data, cmap='magma', norm=norm, origin='lower',
               zorder=-50)  # negative zorder to be below pysiaf aperture fill zorder


    #-- Annotate coordinate grid
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='white', ls='dotted', alpha=0.5)

    #-- Mark guide star
    gscolor='yellow'
    v1color='white'
    Jcolor = 'chartreuse'

    slew = visit.slew
    # Compute expected V3PA
    v3pa_at_gs = visit.slew.GSPA + (0 if visit._no_gspa_yoffset else fgs_aperture.V3IdlYAngle)

    guidemode = slew.GUIDEMODE

    if v1_ra > j_ra:
        v1_ha, j_ha = 'right', 'left'
    else:
        v1_ha, j_ha = 'left', 'right'

    # subsequent annotations can mess up the axes limits, so save here and restore later
    # this is a hacky workaround
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    plt.scatter([v1_ra], [v1_dec], marker='+', s=50, color=v1color,
            transform=ax.get_transform('icrs'))
    plt.text(v1_ra, v1_dec, "\nV1 axis",
             transform=ax.get_transform('icrs'), fontsize=10,
             horizontalalignment=v1_ha, verticalalignment='top', color=v1color)

    plt.scatter([j_ra], [j_dec], marker='*', s=50, color=Jcolor,
            transform=ax.get_transform('icrs'))
    plt.text(j_ra, j_dec, "\nJ axis",
             transform=ax.get_transform('icrs'), fontsize=10,
             horizontalalignment=j_ha, verticalalignment='top', color=Jcolor)



    #-- Plot JWST apertures at that visit's orientation(s)

    if guidemode=='COARSE':
        # we only have a science attitude
        attmatsci = visit.get_attitude_matrix(step='slew')
        guide_det_info = ""
        gslabel = "\n'pseudo guide star'\n(slew coordinates\n for Coarse pointing)"
        plt.scatter(visit.slew.GSRA, visit.slew.GSDEC, s=160, edgecolor=gscolor, facecolor='none',
                    transform=ax.get_transform('icrs'))
        plt.text(visit.slew.GSRA, visit.slew.GSDEC, gslabel,
                 transform=ax.get_transform('icrs'),
                 horizontalalignment='left', verticalalignment='top', color=gscolor)
        plt.text(0.02, 0.07,
             f"Pseudo guide star at {visit.slew.GSRA:.7f}, {visit.slew.GSDEC:.7f}, GSPA={visit.slew.GSPA}",
             color=gscolor, transform=ax.transAxes, verticalalignment='bottom')

    else:
        # There are possibly distinct science and ID attitudes
        # And there may be multiple guide star
        attmatsci = visit.get_attitude_matrix(step='sci')
        attmatid = visit.get_attitude_matrix(step='id')
        # print("ID ATTITUDE:", attmatid)
        # print("SCIENCE ATTITUDE:", attmatsci)
        fgs_detector = 1 if visit.slew.DETECTOR=='GUIDER1' else 2

        # Plot FGS aperture at ID attitude
        fgs_aperture.set_attitude_matrix(attmatid)
        fgs_aperture.plot(frame='sky', transform=ax.get_transform('icrs'), color=gscolor, fill=False)
        plt.text(0.02, 0.02, f"Yellow = ID attitude",
                color=gscolor, transform=ax.transAxes, horizontalalignment='left')

        # Set up some variables we'll use below
        gs_warning_text = ""
        gs_linestyles = ('-', '--', ':', 'dashdot')

        # Plot (and check) the different possible GS reference stars

        # first, check if all guiding is on same guider
        guiders_used = set()
        for gs_id, gs in enumerate(visit.guide_activities):
            guiders_used.add(gs.DETECTOR)
        single_guider = len(guiders_used)==1
        guide_det_info = f", on FGS{fgs_detector}" if single_guider else ", G1 or G2"

        for gs_id, gs in enumerate(visit.guide_activities):
            if gs.args[1]=='FGSVERMAIN': continue


            # For each FGSMAIN call, work out the attitude matrix for the ID attitude:
            attmatid = visit.get_attitude_matrix(step='id', gscandidate=gs_id + 1)
            if verbose:
                print(f"GS {gs_id + 1}, ID attmat: {attmatid}")
            fgs_aperture.set_attitude_matrix(attmatid)

            # how many reference stars does this guide star have?
            nrefs = max([n for n in range(0,11) if hasattr(gs, f'REF{n}X')]+[0])

            # Compute the RA, Dec of the guide star. We do this using ID X,Y coordinates because
            # second and subsequent FGSMAIN calls seem not to always have GSRA, GSDEC, but they
            # always have GSXID, GSYID
            gs_radec = fgs_aperture.idl_to_sky(gs.GSXID, gs.GSYID)
            try:
                gspa = gs.GSPA
            except AttributeError:
                gspa = visit.slews[gs_id].GSPA

            # Plot them, using shades fading out for later attempts
            alpha = max(1-gs_id*0.2, 0.5)
            gs_ls = gs_linestyles[min(gs_id, len(gs_linestyles)-1)]

            fgs_aperture.plot(frame='sky', transform=ax.get_transform('icrs'), color=gscolor, ls=gs_ls, fill=False, alpha=alpha)
            plt.scatter(*gs_radec, s=160, edgecolor=gscolor, facecolor='none',
                        transform=ax.get_transform('icrs'), alpha=alpha)
            plt.text(*gs_radec, f"   \n   GS candidate {gs_id+1}, with {nrefs} ref", alpha=alpha,
                    color=gscolor, transform=ax.get_transform('icrs'), horizontalalignment='left')
            extra_gs_text = "" if single_guider else f" on {gs.DETECTOR.replace('UIDER','')}"
            plt.text(0.02, 0.07 + 0.02*(len(visit.guide_activities) - gs_id - 1),
                     f"Guide star candidate {gs_id+1} at {gs_radec[0]:.7f}, {gs_radec[1]:.7f}, GSPA={gspa}{extra_gs_text}",
                     color=gscolor, alpha=alpha, transform=ax.transAxes, verticalalignment='bottom')

            # Iterate over reference stars for this guide star. Plot them, and sanity check their distinctness
            for ref_id in range(1,10):
                if hasattr(gs, f'REF{ref_id}X'):
                    refx = getattr(gs, f'REF{ref_id}X')
                    refy = getattr(gs, f'REF{ref_id}Y')
                    refradec = fgs_aperture.idl_to_sky(refx, refy)

                    # Overplot
                    plt.scatter(*refradec, marker='D', s=40, color='orange', alpha=alpha, facecolor='none',
                                transform=ax.get_transform('icrs'))
                    if ref_id==1 and gs_id==0:
                        plt.text(*refradec, "GS References:  ",
                                 transform=ax.get_transform('icrs'),
                                 horizontalalignment='right', verticalalignment='center', color='orange')

                    # Let's do some sanity checks
                    # Is this reference star distinct from the guide star?
                    gs_ref_separation = np.sqrt((gs.GSXID-refx)**2 + (gs.GSYID - refy)**2 )
                    if gs_ref_separation < 2:
                        errmsg =f"WARNING: FOR GS #{gs_id+1}, GUIDE STAR AND REF {ref_id} ARE LESS THAN 2 ARCSEC APART\nMAY NOT BE WELL SEPARATED ENOUGH\n"
                        gs_warning_text += errmsg
                    # Is this reference star distinct from the other references?
                    for other_ref_id in range(1,ref_id):
                        other_refx = getattr(gs, f'REF{other_ref_id}X')
                        other_refy = getattr(gs, f'REF{other_ref_id}Y')
                        ref_ref_separation = np.sqrt((other_refx-refx)**2 + (other_refy - refy)**2 )
                        if ref_ref_separation < 2:
                            errmsg = f"WARNING: FOR GS #{gs_id+1}, REF {ref_id} AND REF {other_ref_id} ARE LESS THAN 2 ARCSEC APART\nMAY NOT BE WELL SEPARATED ENOUGH\n"
                            gs_warning_text += errmsg

        if gs_warning_text != "":
            print(gs_warning_text)
            plt.gcf().text(0.5, 0.1, gs_warning_text,
                     horizontalalignment='center',
                     color='red', fontweight='bold', fontsize=18, zorder=1000)
        #plot_gs_id_references(visit.s

    # Plot all apertures, faintly
    pysiaf.siaf.plot_main_apertures(frame='sky', darkbg=True,
                               attitude_matrix=attmatsci, transform=ax.get_transform('icrs'), alpha=0.4, fill=False)

    # Plot the active apertures, more visibly
    for apername in visit.apertures_used():
        # look up aperture from that aperture name

        aper_key = apername[0:4] if apername.startswith("M") else apername[0:3]  # NRC, NRS, NIS, FGS, or MIRI

        aperture = SIAFS[aper_key][apername]

        # plot at the correct attitude
        aperture.set_attitude_matrix(attmatsci)
        aperture.plot(frame='sky', transform=ax.get_transform('icrs'), color='cyan', fill=True, fill_alpha=0.2)

    # Highlight main WFSC aperture, NIRCam A3
    # TODO make this smart about module A vs B? And about other templates
    if 'WFSC' in visit.template:
        nrca3_aperture = pysiaf.Siaf("NIRCam").apertures['NRCA3_FULL']
        nrca3_aperture.set_attitude_matrix(attmatsci)
        nrca3_aperture.plot(frame='sky', transform=ax.get_transform('icrs'), color='white', fill=False)

    # Does this visit have multiple pointings for a mosaic? If so, plot those.
    mosaic_pos = []
    if visit.uses_sams():
        gsoffset = np.zeros(2)  # start with 0,0 offsets initially at start of visit, relative to initial science attitude
        attmatsci_mosaic = visit.get_attitude_matrix(step='sci', fgs_delta_from_sam=gsoffset)
        for iact, act in enumerate(visit.si_activities):
            if act.scriptname == 'SCSAMMAIN':
                if hasattr(act, "FGS1DELTAX"):
                    sam = (act.FGS1DELTAX, act.FGS1DELTAY)
                else:
                    sam = (act.FGS2DELTAX, act.FGS2DELTAY)
                gsoffset[0] += sam[0]
                gsoffset[1] += sam[1]
                mosaic_pos.append(gsoffset.copy())
                if verbose:
                    print(f"Mosaic offset SAM: {sam}\tcumulative: {gsoffset} arcsec")

                attmatsci_mosaic = visit.get_attitude_matrix(step='sci', fgs_delta_from_sam=gsoffset)
            else:
                for apername in visit.apertures_used():
                    if apername.startswith('FGS'): continue # don't plot these for each pointing
                    aper_key = apername[0:4] if apername.startswith("M") else apername[
                                                                              0:3]  # NRC, NRS, NIS, FGS, or MIRI
                    aperture = SIAFS[aper_key][apername]
                    aperture.set_attitude_matrix(attmatsci_mosaic)
                    aperture.plot(frame='sky', transform=ax.get_transform('icrs'), ls='--',
                                  color='cyan', fill_color='cyan', fill=True, fill_alpha=0.02, alpha=0.5)
    else:
        n_mosaic_pos = 0


    plt.text(0.02, 0.98, f"Pointing for {os.path.basename(visit.filename)}", color='white', fontweight='bold',
            transform=ax.transAxes, fontsize=12, verticalalignment='top')
    if len(mosaic_pos) > 0:
        mosaic_pos = np.asarray(mosaic_pos)
        plt.text(0.03, 0.95, f"with {len(mosaic_pos)} mosaic pointings, spanning {np.ptp(mosaic_pos[:,0]):.1f}, {np.ptp(mosaic_pos[:,1]):.1f} arcsec in FGS X,Y", color='cyan', fontweight='bold',
                 transform=ax.transAxes, fontsize=11, verticalalignment='top')
    template_ypos = 0.94 if len(visit.template) > 35 else 0.98   # avoid text overlap for very long names
    plt.text(0.98, template_ypos, f"Template = {visit.template}", color='white',
            transform=ax.transAxes, horizontalalignment='right',
            verticalalignment='top')
    plt.text(0.02, 0.05, f"Guide mode = {guidemode}{guide_det_info}",
            color=gscolor, transform=ax.transAxes)
    plt.text(0.98, 0.05, f"Shaded apertures are used\n in this observation",
            color='cyan', transform=ax.transAxes, horizontalalignment='right')
    plt.text(0.98, 0.02, f"Cyan = Science attitude",
            color='cyan', transform=ax.transAxes, horizontalalignment='right')
    plt.text(0.36, 0.02, f"V1 axis at {v1_ra:.3f}, {v1_dec:.3f}\nV3PA={v3pa_at_gs:.3f} at GS, {v3pa_at_v1:.3f} at V1 axis",
            color='white', transform=ax.transAxes)

    apt_program_id = int(os.path.basename(visit.filename)[1:6])
    if verbose:
        print(f"APT program: {apt_program_id}")
    if apt_program_id in PROGRAMS_WITH_SEGMENT_GUIDING and visit.uses_guiding():
        plt.text(0.02, np.max( [0.125, 0.07 + 0.02*(len(visit.guide_activities) +1)])
                 , f"Segment guiding may be used in this program\nThe guide star indicated may be a segment PSF offset from the star location",
            color=gscolor, transform=ax.transAxes)



    if visit._no_gspa_yoffset:
        plt.text(0.00, -0.035,
                 f"Interpreting GSPA parameter like PPS$\\leq$14.14.1:\nGSPA does not include FGS Yics rotation.",
                 color='orangered', transform=ax.transAxes, verticalalignment='top')

    # re-establish the axes limits
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    return ax


def plot_mosaic_pointings_fov(visitlist, center_on_visit=0, subplotspec=None, verbose=False, use_dss=False, title=None,
                              crop_for_no12=False):
    """Plot multiple visits together in one FOV as a mosaic

    CAUTION: Only tested so far on an example mock OTE-01 visit set from Normal Ops 12,
    which is a special case for NIRCam mosaic in coarse point.
    This is not (yet..?) a generalized mosaic viewer! YMMV.

    """
    # Use one visit, by default the first visit in the list, to set up the plot FOV and image
    visit = visitlist[center_on_visit]

    fgs_aperture_name = visit.get_guider_aperture(return_name=True)
    fgs_aperture = SIAFS['FGS'].apertures[fgs_aperture_name]
    attmat = visit.get_attitude_matrix(step='sci')
    fgs_aperture.set_attitude_matrix(attmat)
    mcf_ra, mcf_dec = fgs_aperture.tel_to_sky(0, -468)  # RA, Dec of master chief ray location (between NIRCam A+B)

    # Retrieve image
    if use_dss:
        img_hdu = retrieve_dss_image(visit, ra=mcf_ra, dec=mcf_dec, verbose=verbose, fov=0.5)
    else:
        img_hdu = retrieve_2mass_image(visit, ra=mcf_ra, dec=mcf_dec, verbose=verbose, fov=0.5)

    wcs = astropy.wcs.WCS(img_hdu[0].header)

    # -- Setup plot and axes
    if subplotspec is None:
        plt.figure(figsize=(16, 9))
        ax = plt.subplot(projection=wcs)
    else:
        ax = plt.subplot(subplotspec, projection=wcs)

    # -- Display the 2MASS image
    norm = astropy.visualization.ImageNormalize(img_hdu[0].data,
                                                interval=astropy.visualization.PercentileInterval(99.99),
                                                stretch=astropy.visualization.AsinhStretch(a=0.001))
    plt.imshow(img_hdu[0].data, cmap='magma', norm=norm, origin='lower',
               zorder=-50)  # negative zorder to be below pysiaf aperture fill zorder

    # use colormap to shade each visit distinctly:
    cmap = matplotlib.cm.cool

    for i, visit in enumerate(visitlist):

        slew = visit.slew
        # Compute expected V3PA
        v3pa_at_gs = visit.slew.GSPA + (0 if visit._no_gspa_yoffset else fgs_aperture.V3IdlYAngle)

        if slew.GUIDEMODE != 'COARSE':
            raise RuntimeError("We only expect coarse point for OTE-01 mosaic tiles")

        # we only have a science attitude
        attmatsci = visit.get_attitude_matrix(step='slew')

        gsoffset = np.zeros(2)  # start with 0,0 offsets initially at start of visit

        centers = []

        for iact, act in enumerate(visit.activities):

            if act.scriptname == 'NRCMAIN':
                # Draw NIRCam apertures
                for apername in visit.apertures_used():

                    if apername.startswith('NRS'):
                        continue  # ignore parallel nirspec darks from NO-12

                    # look up aperture from that aperture name

                    aper_key = apername[0:4] if apername.startswith("M") else apername[
                                                                              0:3]  # NRC, NRS, NIS, FGS, or MIRI

                    aperture = SIAFS[aper_key][apername]

                    col = cmap(i / len(visitlist))
                    # plot at the correct attitude
                    aperture.set_attitude_matrix(attmatsci)
                    centers.append(aperture.sci_to_sky(1024, 1024))

                    aperture.plot(frame='sky', transform=ax.get_transform('icrs'),
                                  color=col, fill_color=col, fill=True, fill_alpha=0.2)

                    if aperture.AperName == 'NRCA3_FULL':
                        c0, c1 = aperture.det_to_sky(1024, 1024)
                        plt.text(c0, c1, f"v{i:03d}a{iact:1d}",
                                 color='white', transform=ax.get_transform('icrs'), horizontalalignment='left')
            elif act.scriptname == 'SCSAMMAIN':
                gsoffset[0] += act.FGS1DELTAX
                gsoffset[1] += act.FGS1DELTAY
                if verbose:
                    print(gsoffset)

                attmatsci = visit.get_attitude_matrix(step='slew', fgs_delta_from_sam=gsoffset)

    if title:
        plt.title(f"Visit file pointings for {title}")

    if crop_for_no12:
        # Hack: hard-coded zoom for NO-12 OTE-01 test
        ax.set_xlim(424, 1024)
        ax.set_ylim(150, 750)

    plt.text(0.01, 0.01, "Labels show visit:activity in each NRCA3 pointing", color='white',
             transform=ax.transAxes)


def plot_gs_id_references(activity_statement):
    """ Plot reference stars used for guide star ID
    To do so we have to convert from detector coordinates to sky coordinates
    """


##--- Functions for plotting the visit field of REGARD

def plot_circle_coords(x1, x2, from_frame, to_frame, ax, **plot_kwargs):
    """Plot a circle (great or otherwise) with astronomical coordinate transforms"""
    # make a circle of coords
    from_plane = coords.SkyCoord(x1, x2, unit='deg', frame=from_frame)
    return plot_coords_in_frame(from_plane, ax, to_frame, **plot_kwargs)

def plot_coords_in_frame(skycoords, ax, to_frame, label=None, labelcolor='black', scatter=False, **plot_kwargs):
    """ Plot in ICRS or Barycentric coordinates, with appropriate transformations.

    Note this implicitly/automatically take care of flipping the sign of the X axis
    to implement the astronomical convention of R.A. increasing to the left, which
    is awkward to do in matplotlib map projections.
    """
    coords_in_frame = skycoords.transform_to(to_frame)
    # pull out what we want to plot
    if to_frame == 'icrs' or to_frame == 'gcrs':
        plot_x = coords_in_frame.ra.wrap_at('180d').radian
        plot_y = coords_in_frame.dec.radian
    elif to_frame == 'geocentricmeanecliptic' or to_frame=='barycentricmeanecliptic':
        plot_x = coords_in_frame.lon.wrap_at('180d').radian
        plot_y = coords_in_frame.lat.radian
    else:
        raise ValueError("Unsupported Frame")

    # note we MUST FLIP THE SIGN FOR X, since matplotlib map projections
    # don't let us make R.A. increase to the left
    plot_function = ax.scatter if scatter else ax.plot
    plot_function(-plot_x, plot_y, **plot_kwargs)

    if label:
        ax.text(-plot_x, plot_y, "\n"+label, color=labelcolor, horizontalalignment='center', verticalalignment='top', zorder=10)
        #print(-plot_x, plot_y, label)



def show_field_of_regard_ecliptic(visit, datetime=None, projection='lambert', subplotspec=None, show_sun=False):
    """Plot JWST field of regard

    Lotsa complicated stuff here with map projections and coordinates...
    """

    if datetime is None:
        # What date/time are we generating this plot for?
        datetime = visit.time_early # astropy.time.Time.now()


    # Determine where the sun will be then.
    # NOTE - this is the apparent position as seen from Earth, and may be off by up to ~0.3 deg from the
    # perspective of JWST in its large orbit around L2. This is not a high accuracy calculation.
    sun = coords.get_sun(datetime).transform_to('geocentricmeanecliptic')
    center_longitude = sun.lon.radian+np.pi
    subplot_kwargs = {'projection': projection}
    if projection=='lambert':
        # Need to flip sign on center longitude, because of R.A. convention
        subplot_kwargs['center_longitude']= -center_longitude
    if subplotspec is None:
        ax = plt.subplot(**subplot_kwargs)
    else:
        ax = plt.subplot(subplotspec, **subplot_kwargs)

    ax.set_xticklabels([])
    ax.grid(True)

    #
    # if 0:
    #     # Now do a fill-between for those two circles.
    #     # This is tedious since we need to convert to Cartesian for fill_between
    #     n=100
    #     x = np.linspace(-ro, ro, n)
    #     y2 = np.sqrt(ro**2-x**2)
    #     y1 = np.zeros_like(x)
    #     w = np.abs(x) < ri
    #     y1[w] = np.sqrt(ri**2-x[w]**2)
    #
    #     plt.fill_between(x+0.5, y1+0.5, y2+0.5, color='green', alpha=0.1, transform=ax.transAxes,)
    #     plt.fill_between(x+0.5, -y1+0.5, -y2+0.5, color='green', alpha=0.1, transform=ax.transAxes,)
    #
    #     plot_circle_coords(np.arange(361), np.zeros(361), 'galactic', 'barycentricmeanecliptic', ax,
    #                        ls='none',marker='.', markersize=2, color='maroon', alpha=0.3)
    #
    #     plot_circle_coords(np.arange(361)-180, np.zeros(361), 'centricmeanecliptic', 'geocentricmeanecliptic', ax,
    #                        ls='none',marker='.', markersize=2, color='black', alpha=0.3)
    #
    #     #ecliptic_plane = coords.SkyCoord(  np.arange(361)-180, np.zeros(361), unit='deg', frame='barycentricmeanecliptic')
    #     #ax.plot(ecliptic_plane.lon.wrap_at('180d').radian, ecliptic_plane.lat.radian, ls='none',
    #     #        marker='.', markersize=2, color='black', alpha=0.3)
    #
    #     # Show the target pointing!
    #     gs = coords.SkyCoord(visit.slew.GSRA, visit.slew.GSDEC, unit='deg')
    #     plot_coords_in_frame(gs, ax, 'barycentricmeanecliptic',  marker='*',color='red', markersize=20,)
    #
    #
    #     plt.plot(antisun.lon.radian, antisun.lat.radian, marker='+', color='black', markersize=10, markeredgewidth=3,)
    # Annotate
    if projection == 'lambert':
        plt.text(0.5, 0.55, 'anti-sun', transform=ax.transAxes, color='black', horizontalalignment='center')
        plt.text(0.5, 0.5+0.38, 'N CVZ', transform=ax.transAxes, color='blue')
        plt.text(0.5, 0.5-0.38, 'S CVZ', transform=ax.transAxes, color='blue', verticalalignment='top')
        plt.text(0.5-0.35, 0.54, 'ram', transform=ax.transAxes, color='dimgray', verticalalignment='bottom')
        plt.text(0.5+0.35, 0.54, 'wake', transform=ax.transAxes, color='dimgray', verticalalignment='bottom', horizontalalignment='right')

    #plt.text(0.99, 0.01, f"JWST field of regard\non {datetime.to_value('iso',subfmt='date')}\n[Ecliptic coords]",
    #         transform=ax.transAxes, horizontalalignment='right')

    plt.title(f"JWST field of regard\non {datetime.to_value('iso',subfmt='date')}\n[Ecliptic coords]",)

    # Plot all the markers
    plot_celestial_markers(visit, ax, 'geocentricmeanecliptic', show_sun=show_sun, datetime=datetime)

    return ax




def show_field_of_regard_ra_dec(visit, datetime=None, subplotspec=None, labelaxes=False, label_markers=False):
    """Plot celestial sphere in regular RA, Dec

    """

    # Setup a map in Mollweide projection and label for RA and Dec.
    if subplotspec is None:
        ax = plt.gcf().add_subplot( projection='mollweide')
    else:
        ax = plt.subplot(subplotspec, projection='mollweide')
    plt.title(f"R.A., Decl.\n[ICRS Equatorial coords]\n")
    if labelaxes:
        plt.ylabel("Declination")
        plt.xlabel("Right Ascension")
    # Replace the standard x ticks (longitude) with R.A., and include sign flip so R.A. increases to left
    plt.xticks(ticks=np.radians([-120, -60, 0, 60, 120, 180]),
               labels=['8$^h$', '4$^h$', '0$^h$', '20$^h$', '16$^h$', '12$^h$'])
    ax.grid(True)

    # Plot all the markers
    # Use GCRS coordinates, which is like ICRS but references to Earth; this only matters substantially for
    # the sun coordinates in this case.
    plot_celestial_markers(visit, ax, 'gcrs', datetime=datetime, label_markers=label_markers)

    return ax


def plot_celestial_markers(visit, ax, frame='gcrs', datetime=None, show_sun=True, label_markers=False, npts=360+1):
    """Main routine to plot celestial markers for sun, target, and various circles

    Parameters
    ----------
    npts : int
        Number of points to rasterize circles into. Increase this if the plot output is not smooth for larger plot sizes.
    """

    if datetime is None and visit is not None:
        # What date/time are we generating this plot for?
        datetime = visit.time_early # astropy.time.Time.now()

    # Determine where the sun will be then.
    # NOTE - this is the apparent position as seen from Earth, and may be off by up to ~0.3 deg from the
    # perspective of JWST in its large orbit around L2. This is not a high accuracy calculation.
    sun = coords.get_sun(datetime).transform_to(frame)

    # Draw the sun, and anti-sun point
    antisun = sun.directional_offset_by(0, 180*u.deg)
    if show_sun: # it's less useful to show this in the lambert ecliptic proj
        plot_coords_in_frame(sun, ax, frame,  marker='o', markersize=20, color='orange', markeredgecolor='orange', zorder=10,
                             label="Sun" if label_markers else None)
    plot_coords_in_frame(antisun, ax, frame,  marker='+', markersize=10, color='black', markeredgewidth=3,
                         label='anti-Sun' if label_markers else None)

    # Draw galactic and ecliptic planes
    # JWST's CVZs are defined by the Barycentric version of the ecliptic plane.
    plot_circle_coords(np.linspace(0,360,npts), np.zeros(npts), 'galactic', frame, ax,
                       ls='none',marker='.', markersize=2, color='maroon', alpha=0.3)
    plot_circle_coords(np.linspace(0,360,npts)-180, np.zeros(npts), 'barycentricmeanecliptic', frame, ax,
                       ls='none',marker='.', markersize=2, color='black', alpha=0.3)
    # Draw the CVZs
    for ecliptic_lat in [85, -85]:
        plot_circle_coords(np.arange(361) - 180, np.zeros(361)+ecliptic_lat, 'barycentricmeanecliptic', frame, ax,
                           ls='none', marker='.', markersize=0.5, color='blue', alpha=0.3)


    # Shade the allowed field of regard.
    # Given hassles with plotting filled shapes in map projections, a relatively straightforward way to do this
    # ends up being to rasterize a grid over the sky, and scatter plot those points densely
    #lat, lon = np.meshgrid(np.linspace(-np.pi, np.pi, 2*npts), np.linspace(-np.pi / 2, np.pi / 2, npts))
    #skymesh = astropy.coordinates.SkyCoord(lat, lon, frame=frame, unit='rad')

    pas = np.linspace(0,360, npts) * u.deg
    for offset in [85, 135]:
        for_edge = sun.directional_offset_by(pas, offset*u.deg)
        plot_coords_in_frame(for_edge, ax, frame,  marker='.',color='green', markersize=1, ls='none')

    for offset in np.linspace(85, 135, npts//3):
        for_circ = sun.directional_offset_by(pas, offset * u.deg)
        plot_coords_in_frame(for_circ, ax, frame, marker='o', color='#E5F2E5', markersize=1, ls='none', zorder=-30)

    # Show the target pointing!
    if visit is not None:
        gs = coords.SkyCoord(visit.slew.GSRA, visit.slew.GSDEC, unit='deg', frame='icrs')
        plot_coords_in_frame(gs, ax, frame,  marker='*',color='red', markersize=20,
                             label='Target' if label_markers else None)

    gal_center = astropy.coordinates.SkyCoord(0, 0, unit='deg', frame='galactic').transform_to('icrs')
    plot_coords_in_frame(gal_center, ax, frame, marker='o', markersize=5, color='maroon',
                         label='galactic center' if label_markers else None, labelcolor='maroon',)

    # Extra markers to plot while debugging transforms
    #origin = coords.SkyCoord(0, 0, unit='deg', frame='icrs')
    #one = coords.SkyCoord(15, 0, unit='deg', frame='icrs')
    #plot_coords_in_frame(origin, ax, frame,  marker='^',color='green', markersize=20,)
    #plot_coords_in_frame(one, ax, frame,  marker='^',color='lightgreen', markersize=20,)
    #two = coords.SkyCoord(-30, 0, unit='deg', frame='icrs')
    #plot_coords_in_frame(two, ax, frame,  marker='^',color='darkgreen', markersize=20,)

def show_pitch_roll(visit, subplotspec_pitch=None, subplotspec_roll=None):
    """ Plot pitch and roll relative to the sun

    We here distinguish between sun pitch and roll, and vehicle pitch and roll.
    These are closely related but not quite the same. See John Isaac's document
    "Conversion between Vehicle Roll and Sun Roll".  Briefly::

        The Sun roll angle S_r is the rotation angle about the vector normal to the Sun vector,
        in the plane defined by the +V1 axis and the Sun vector, and in the -V1 axis direction.

        The vehicle roll angle V_R is the rotation about the +V1 axis from the optimal vehicle orientation.

        The vehicle sun angle A_s is the angle from the +V1 axis to the sun (i.e. typically ranging from 85
        to 135 degrees). The vehicle pitch angle is then Vp = 90 - A_s (i.e. typically from -5 to +45)

        The Sun pitch angle S_p is the angle from the -V3 axis to the projection of the Sun vector in the V1-V3 plane,
        measured with respect to the +V2 axis (i.e. typically ranging from +5 to -45 degrees).

    """
    import jwst_gtvt.jwst_tvt
    #import jwst_gtvt.find_tgt_info

    if subplotspec_pitch is None or subplotspec_roll is None:
        fig, subplots = plt.subplots(nrows=2, figsize=(16,9))
        ax_pitch, ax_roll = subplots
    else:
        ax_pitch = plt.subplot(subplotspec_pitch)
        ax_roll = plt.subplot(subplotspec_roll)

    dt = visit.time_late - visit.time_early
    nt = 30

    times = visit.time_early + dt * np.linspace(0, 1, nt)
    gs = astropy.coordinates.SkyCoord(visit.slew.GSRA, visit.slew.GSDEC, unit='deg',
                                      frame='icrs')


    # find where the V1 axis will be for this visit
    attmat = visit.get_attitude_matrix(step='slew')
    fgs_aperture_name = visit.get_guider_aperture(return_name=True)
    fgs_aperture = SIAFS['FGS'].apertures[fgs_aperture_name]
    fgs_aperture.set_attitude_matrix(attmat)
    #v1coords = astropy.coordinates.SkyCoord(*fgs_aperture.tel_to_sky(0,0), unit='deg', frame='icrs')
    v1ra, v1dec = fgs_aperture.tel_to_sky(0,0)
    v1ra = np.deg2rad(v1ra)
    v1dec = np.deg2rad(v1dec)


    #sun = astropy.coordinates.get_sun(times)

    # Load JWST ephemeris from jwst_gtvt
    A_eph = jwst_gtvt.jwst_tvt.Ephemeris()

    v3pa_nominal = np.zeros(nt, float)
    max_allowed_roll = np.zeros(nt, float)
    sun_angle = np.zeros(nt, float)


    for i in range(nt):
        # the following all in radians, for compatibility with jwst_gtvt internals
        atime = times[i]
        try:
            v3pa_nominal[i] = A_eph.normal_pa(atime.mjd, v1ra, v1dec)
            # Determine apparent sun position at that time.
            sun_ra, sun_dec = A_eph.sun_pos(atime.mjd)

            sun_angle[i] = jwst_gtvt.find_tgt_info.angular_sep(sun_ra, sun_dec, v1ra, v1dec)

            max_allowed_roll[i] = jwst_gtvt.find_tgt_info.allowed_max_vehicle_roll(sun_ra, sun_dec, v1ra, v1dec)
        except:
            sun_angle[i] = np.nan
            max_allowed_roll[i] = np.nan

    # convert to degrees
    sun_angle = np.rad2deg((sun_angle))
    max_allowed_roll = np.rad2deg((max_allowed_roll))
    v3pa_nominal = np.rad2deg((v3pa_nominal))

    vehicle_pitch = sun_angle - 90

    V3PA = visit.slew.GSPA + fgs_aperture.V3IdlYAngle
    vehicle_roll = V3PA - v3pa_nominal   # TODO update with V3PA not GSPA

    vr = np.deg2rad(vehicle_roll)
    vp = np.deg2rad(vehicle_pitch)
    sun_roll = np.rad2deg(np.arcsin( np.sin(vr) * np.cos(vp) ))
    sun_pitch = np.rad2deg(np.arctan( np.tan(vp) / np.cos(vr)))

    #print( vehicle_pitch.mean(), sun_pitch.mean(), vehicle_roll.mean(), sun_roll.mean())

    with astropy.visualization.time_support():

        ax_pitch.plot_date(times.plot_date, vehicle_pitch, ls='-', lw=3, marker=None, label='Vehicle pitch')
        ax_pitch.plot_date(times.plot_date, sun_pitch, ls='--', lw=1, color='orange', marker=None, label='Sun pitch')

        ax_pitch.set_ylim(-10, 50)
        for lim in [-5, 45]:
            ax_pitch.axhline(lim, ls='--', color='red')
        ax_pitch.set_ylabel("Pitch\n[deg]")
        ax_pitch.set_yticks([-5, 0, 22.5, 45])
        ax_pitch.axhline(0, color='black', alpha=0.3, ls=":")
        # ax_pitch.set_xlabel("Date time [UTC]")

        ax_roll.plot_date(times.plot_date, vehicle_roll, ls='-', lw=3, marker=None, label='Vehicle roll')
        ax_roll.plot_date(times.plot_date, sun_roll, ls='--', lw=1, color='orange', marker=None, label='Sun roll')
        ax_roll.set_ylim(-8, 8)

        ax_roll.plot_date(times.plot_date, max_allowed_roll, ls='--', color='red', marker=None)
        ax_roll.plot_date(times.plot_date, -max_allowed_roll, ls='--', color='red', marker=None)

        ax_roll.set_ylabel("Roll\n[deg]")
        ax_roll.set_xlabel("Date time [UTC]")

        for ax in [ax_pitch, ax_roll]:
            for label in ax.get_xticklabels():
                label.set_rotation(40)
            ax.legend(loc='upper right', fontsize=8)

    ax_pitch.set_title(f"Mean pitch: {vehicle_pitch.mean():.2f}$^\\circ$, Roll range: {vehicle_roll[0]:.2f}$^\\circ$ to {vehicle_roll[-1]:.2f}$^\\circ$",
                       fontweight='bold', fontsize=9, color='C0')

    if np.any( np.abs(vehicle_roll) > max_allowed_roll):
        plt.text(0.98, 0.1, 'WARNING! Vehicle roll may be outside of limit at some times?\n(IGNORE this warning; need to update with flight ephemeris)', color='red',
                 #fontweight='bold',
                 horizontalalignment='right',fontsize=15,
                 transform=plt.gcf().transFigure)


# Functions for combining multiple plot panels into one page:


def multi_plot(visit, verbose=False, save=False, use_dss=False, no_gspa_yoffset=False, output_dir=None, center_visit=None):
    """ Main top-level function for visitviewer"""
    fig = plt.figure(figsize=(16, 9))

    # Set up a bunch of plot axes, complicatedly via nested grids
    # (this could be simplified but I'm adapting from other existing code...)
    gs_outer = plt.GridSpec(1, 2, width_ratios=[9,7])
    gs_outer.update(bottom=0.05, top=0.95, wspace=0.2, left=0.04, right=0.98)

    r_gridspec_kw = {'hspace': 0.2, 'wspace': 0.2, 'height_ratios': [0.3, 1, 1, 0.3], 'width_ratios': [1, 1.5]}
    gs_r = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=gs_outer[1],
                                            **r_gridspec_kw)
    r2_gridspec_kw = {'hspace': 0.05, 'wspace': 0.2, 'height_ratios': [0.3, .7, 0.3, 0.3, 0.4]}
    gs_r2 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs_outer[1],
                                             **r2_gridspec_kw)


    # handle edge case: visit with no slew
    if visit.slew._is_no_slew:
        plt.text(0.5, 0.5, f"Unpointed visit! No slew defined\n(OSS SCNOSLEWMAIN)\nCannot plot coordinates or field of view meaningfully, sorry!",
             color='red',
             fontsize=15, horizontalalignment='center', fontweight='bold',
             transform=fig.transFigure)


    else:
        # Now make the plots, via the functions defined above
        plot_visit_fov(visit, subplotspec=gs_outer[0], use_dss=use_dss, verbose=verbose, center_visit=center_visit)
        show_field_of_regard_ecliptic(visit, subplotspec=gs_r[1, 0])
        show_field_of_regard_ra_dec(visit, subplotspec=gs_r[1, 1])
        show_pitch_roll(visit, gs_r2[2, 0], gs_r2[3, 0])

    # Annotate labels
    plt.text(0.9, 0.05, f"Visit times EARLY =  {visit.time_early} UTC\n"
                        f"            LATE =   {visit.time_late} UTC",
             #                    f"            CUTOFF = {visit.time_cutoff}",
             fontsize=12, horizontalalignment='right', fontweight='bold',
             transform=fig.transFigure)

    plt.text(0.55, 0.95, visit.short_summary(),
             fontsize = 12, horizontalalignment = 'left', fontweight = 'bold',
             verticalalignment='top',
             transform = fig.transFigure)

    if save:
        outname = f"{visit.visitid}_view.pdf"
        if output_dir is not None:
            outname = os.path.join(output_dir, outname)
        plt.savefig(outname)
        print(f"File saved to {outname}")

        if platform.system()=='Darwin':
            subprocess.run(["open", outname])



def mosaic_plot(visitlist, verbose=False, save=False, use_dss=False, no_gspa_yoffset=False, output_dir=None):
    """ Main top-level function for visitviewer, special version for mosaics"""
    fig = plt.figure(figsize=(16, 9))

    # Set up a bunch of plot axes, complicatedly via nested grids
    # (this could be simplified but I'm adapting from other existing code...)
    gs_outer = plt.GridSpec(1, 2, width_ratios=[9,7])
    gs_outer.update(bottom=0.05, top=0.95, wspace=0.2, left=0.04, right=0.98)

    r_gridspec_kw = {'hspace': 0.2, 'wspace': 0.2, 'height_ratios': [0.3, 1, 1, 0.3], 'width_ratios': [1, 1.5]}
    gs_r = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=gs_outer[1],
                                            **r_gridspec_kw)
    r2_gridspec_kw = {'hspace': 0.05, 'wspace': 0.2, 'height_ratios': [0.3, .7, 0.3, 0.3, 0.4]}
    gs_r2 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs_outer[1],
                                             **r2_gridspec_kw)

    # Use the first visit as reference for field of regard and pitch/roll
    visit=visitlist[0]

    # Now make the plots, via the functions defined above
    plot_mosaic_pointings_fov(visitlist, subplotspec=gs_outer[0], use_dss=use_dss, verbose=verbose)
    show_field_of_regard_ecliptic(visit, subplotspec=gs_r[1, 0])
    show_field_of_regard_ra_dec(visit, subplotspec=gs_r[1, 1])
    show_pitch_roll(visit, gs_r2[2, 0], gs_r2[3, 0])

    # Annotate labels
    plt.text(0.9, 0.05, f"Visit times EARLY =  {visit.time_early} UTC\n"
                        f"            LATE =   {visit.time_late} UTC",
             #                    f"            CUTOFF = {visit.time_cutoff}",
             fontsize=12, horizontalalignment='right', fontweight='bold',
             transform=fig.transFigure)

    plt.text(0.55, 0.95, visit.short_summary(),
             fontsize = 12, horizontalalignment = 'left', fontweight = 'bold',
             verticalalignment='top',
             transform = fig.transFigure)

    if save:
        outname = f"{visit.visitid}_mosaic_view.pdf"
        if output_dir is not None:
            outname = os.path.join(output_dir, outname)
        plt.savefig(outname)
        print(f"File saved to {outname}")

        if platform.system()=='Darwin':
            subprocess.run(["open", outname])
        return outname



def plot_field_of_regard_on_date(date,
                                 use_ecliptic_coords=False,
                                 mark_sun_pitch=None,
                                 label_circles=True, label_axes = True,
                                 label_ram_wake=True,
                                 mark_maz=False,
                                 figsize=(12, 6.75),
                                 save=True):
    """Plot JWST field of regard on a specified date.

    Creates an all-sky plot in Mollweide projection, and annotates it with the
    sun, several major celestial circles, and the JWST field of regard on that date.

    Useful for stand-alone plots without any specific visit file yet.

    Parameters
    ----------
    date : str or astropy.time.Time
        Date, as "YYYY-MM-DD", or some other format convertable to an astropy Time object
    use_ecliptic_coords : bool
        Default is to plot in celestial coordinates (R.A., Dec). If this is set,
        show the plot in ecliptic coordinates instead.
    mark_sun_pitch : float
        If set, draw a circle showing the annulus of the specified sun pitch.
        e.g. mark_sun_pitch=-20 shows the peak power attitude locations.
        Note this uses the ACS sign convention in which the field of regard goes from +5 to -45.
    label_circles : bool
        Add text labels for ecliptic, galactic, and cvz marker circles
    label_axes : bool
        Add text labels for RA, Dec axes
    label_ram_wake : bool
        Add text labels for ram and wake directions
    figsize : tuple of floats
        Figure size for Matplotlib

    """

    plt.figure(figsize=figsize)
    npts = 360*4+1  # Higher than the default, for nice smooth plot output

    datetime = astropy.time.Time(date)

    if use_ecliptic_coords:
        frame = 'geocentricmeanecliptic'
        frame_name, xlabel, ylabel = 'Ecliptic', "Ecliptic Longitude", "Ecliptic Latitude"
        yticklabels = [f'{d}$^\\circ$' for d in [120, 60, 0, 300, 240, 180]]
    else:
        # Use GCRS coordinates, which is like ICRS but references to Earth; this only matters substantially for
        # the sun coordinates in this case.
        frame = 'gcrs'
        frame_name, xlabel, ylabel = 'ICRS Equatorial', "Right Ascension", "Declination"
        yticklabels = ['8$^h$', '4$^h$', '0$^h$', '20$^h$', '16$^h$', '12$^h$']

    # Setup a map in Mollweide projection and label for RA and Dec.
    ax = plt.gcf().add_subplot( projection='mollweide')
    ax.set_title(f"All Sky\n[{frame_name} coords]\n")
    if label_axes:
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)

    # Replace the standard x ticks (longitude) with R.A., and include sign flip so R.A. increases to left

    plt.xticks(ticks=np.radians([-120, -60, 0, 60, 120, 180]),
               labels=yticklabels)
    ax.grid(True)


    # Plot all the markers for sun, ecliptic plane, galactic plane
    plot_celestial_markers(visit=None, ax=ax, frame=frame, datetime=datetime, label_markers=True, npts=npts)

    # Optional, mark a circle at some specified sun pitch
    sun = coords.get_sun(datetime).transform_to(frame)
    if isinstance(mark_sun_pitch, bool):
        raise ValueError("The mark_sun_pitch parameter should be a numeric value between 5 and -45 degrees, not a boolean True/False.")
    if mark_sun_pitch is not None:
        pas = np.linspace(0,360, 721) * u.deg
        for offset in [90 - mark_sun_pitch]:
            pitch_circle = sun.directional_offset_by(pas, offset*u.deg)
            plot_coords_in_frame(pitch_circle, ax, frame,  marker='.',color='magenta', markersize=1, ls='none')
        # Call again for just a single point to add a label.
        plot_coords_in_frame(pitch_circle[200], ax, frame,  marker='.',color='magenta', markersize=1,
                             ls='none', label=f'Sun Pitch {mark_sun_pitch}', labelcolor='magenta')

    # Optional, add additional text labels for features
    if label_circles:
        for cvz_sign in [1, -1]:
            cvz = coords.SkyCoord(0, 90*cvz_sign, unit=u.deg, frame='barycentricmeanecliptic')
            plot_coords_in_frame(cvz, ax, frame, marker='none', label='CVZ', labelcolor='blue')
        # TODO this next bit is a bit awkward, with manual label positioning per coordinate frame
        if use_ecliptic_coords:
            ax.text(np.pi/8, 0.02,'Ecliptic', rotation=0)
            ax.text(np.pi/4, np.pi*0.18,'Galactic plane', color='maroon', rotation=-35)
        else:
            ax.text(np.pi/8, -np.pi/15,'Ecliptic', rotation=-22)
            ax.text(np.pi/6, np.pi*0.2,'Galactic plane', color='maroon', rotation=-40)


    # Mark the ram and wake directions.
    if label_ram_wake:
        sun_ecliptic = sun.transform_to('geocentricmeanecliptic')
        for ramwake_label, ramwake_offset in (('Ram', -90), ('Wake', 90)):
            ramwake_dir = coords.SkyCoord(sun_ecliptic.lon+ramwake_offset*u.deg, sun_ecliptic.lat,
                                          frame='geocentricmeanecliptic').transform_to(frame)
            plot_coords_in_frame(ramwake_dir, ax, frame,  marker='none',
                                 markersize=20, color='orange', markeredgecolor='darkorange', zorder=10,
                                 label=ramwake_label + "\ndirection", labelcolor='darkorange')
            if ramwake_label == 'Ram' and mark_maz:
                # Plot MAZ boundary
                pas = np.linspace(0,360, 721) * u.deg
                with warnings.catch_warnings():
                    # ignore an astropy NonRotationTransformationWarning, which comes up for some reason
                    warnings.simplefilter('ignore')
                    maz_circle = ramwake_dir.transform_to(frame).directional_offset_by(pas, 60*u.deg)
                    maz_dist_sun = sun_ecliptic.separation(maz_circle)
                    # limit to the part of the MAZ circle boundary that is within the FOR
                    maz_circle_in_for = maz_circle[(maz_dist_sun>=85*u.deg) & (maz_dist_sun<135*u.deg)]
                    plot_coords_in_frame(maz_circle_in_for, ax, frame,  marker='.',
                                        markersize=2,  color='darkorange', zorder=1, linestyle='none')

                    # some jumping through hoops to find a nice place to put the text label
                    better_side = np.abs(maz_circle.dec) < 45*u.deg
                    maz_circ_midpoint = np.argmin(np.abs(maz_dist_sun[better_side] - 110*u.deg))
                    plot_coords_in_frame(maz_circle[better_side][maz_circ_midpoint], ax, frame,  marker='none',
                                     markersize=1, color='orange', markeredgecolor='darkorange', zorder=10,
                                     label="MAZ\nboundary", labelcolor='darkorange')




    plt.suptitle(f"JWST Field of Regard on {datetime.iso[0:10]}")
    plt.tight_layout()
    if save:
        fn = f'jwst_field_of_regard_{datetime.iso[0:10]}{"_ecliptic" if use_ecliptic_coords else ""}.png'
        plt.savefig(fn, dpi=150)
        print(f"Plot saved to {fn}")
    return ax
