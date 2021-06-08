import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from urllib.parse import quote
import astropy, astropy.wcs, astropy.visualization
import astropy.io.fits as fits
import astropy.coordinates as coords
import astropy.units as u
import pysiaf
import jwst_gtvt.find_tgt_info
import jwst_gtvt.ephemeris_old2x as EPH




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



##--- Functions for plotting the visit field of VIEW

def retrieve_2mass_image(visit, ra=None, dec=None, verbose=True, redownload=False, filter='K'):
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
    fov = 0.35

    visitname = os.path.splitext(os.path.basename(visit.filename))[0]

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


def plot_visit_fov(visit, verbose=False, subplotspec=None):
    """Make a nice annotated plot of a visit FOV"""

    # let's center the plot on the master chief ray (MCF; between NIRCams).
    # This works better than centering on the guide star.
    # Use SIAF transforms to infer the MCF pointing in RA, Dec.
    fgs_aperture = visit.get_guider_aperture()
    attmat = visit.get_attitude_matrix(step='slew')
    fgs_aperture.set_attitude_matrix(attmat)
    mcf_ra, mcf_dec = fgs_aperture.tel_to_sky(0,-468) # master chief ray location relative to V2V3

    img_hdu = retrieve_2mass_image(visit, ra = mcf_ra, dec=mcf_dec, verbose=verbose)

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

    slew = visit.slew

    plt.scatter(slew.GSRA, slew.GSDEC,  s=200, edgecolor=gscolor, facecolor='none',
            transform=ax.get_transform('icrs'))

    plt.text(slew.GSRA, slew.GSDEC,"\nguide star",
             transform=ax.get_transform('icrs'),
             horizontalalignment='left', verticalalignment='top', color=gscolor)

    # subsequent annotations can mess up the axes limits, so save here and restore later
    # this is a hacky workaround
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()


    #-- Plot JWST apertures at that visit's orientation(s)
    guidemode = slew.GUIDEMODE

    if guidemode=='COARSE':
        # we only have a science attitude
        attmatsci = visit.get_attitude_matrix(step='slew')
        guide_det_info = ""
    else:
        # there are possibly distinct science and ID attitudes
        attmatsci = visit.get_attitude_matrix(step='sci')
        attmatid = visit.get_attitude_matrix(step='id')
        fgs_detector = 1 if visit.slew.DETECTOR=='GUIDER1' else 2
        guide_det_info = f", on FGS{fgs_detector}"

    # Plot FGS aperture at ID attitude (if in a mode using FGS)
    if guidemode != 'COARSE':
        fgs_aperture.set_attitude_matrix(attmatid)
        fgs_aperture.plot(frame='sky', transform=ax.get_transform('icrs'), color=gscolor, fill=False)
        plt.text(0.02, 0.02, f"Yellow = ID attitude",
                color=gscolor, transform=ax.transAxes, horizontalalignment='left')

    # Plot all apertures, faintly
    pysiaf.siaf.plot_main_apertures(frame='sky', darkbg=True,
                               attitude_matrix=attmatsci, transform=ax.get_transform('icrs'), alpha=0.4, fill=False)

    # Plot the active apertures, more visibly
    siafs = load_all_siafs()
    for apername in visit.apertures_used():
        # look up aperture from that aperture name
        aperture = siafs[apername[0:3]][apername]
        # plot at the correct attitude
        aperture.set_attitude_matrix(attmatsci)
        aperture.plot(frame='sky', transform=ax.get_transform('icrs'), color='cyan', fill=True, fill_alpha=0.2)

    # Highlight main WFSC aperture, NIRCam A3
    # TODO make this smart about module A vs B? And about other templates
    if 'WFSC' in visit.template:
        nrca3_aperture = pysiaf.Siaf("NIRCam").apertures['NRCA3_FULL']
        nrca3_aperture.set_attitude_matrix(attmatsci)
        nrca3_aperture.plot(frame='sky', transform=ax.get_transform('icrs'), color='white', fill=False)


    plt.text(0.02, 0.98, f"Pointing for {os.path.basename(visit.filename)}", color='white', fontweight='bold',
            transform=ax.transAxes, fontsize=12, verticalalignment='top')
    template_ypos = 0.94 if len(visit.template) > 35 else 0.98   # avoid text overlap for very long names
    plt.text(0.98, template_ypos, f"Template = {visit.template}", color='white',
            transform=ax.transAxes, horizontalalignment='right',
            verticalalignment='top')
    plt.text(0.02, 0.05, f"Guide star at {visit.slew.GSRA}, {visit.slew.GSDEC}, GSPA={visit.slew.GSPA}\nGuide mode = {guidemode}{guide_det_info}",
            color=gscolor, transform=ax.transAxes)
    plt.text(0.98, 0.05, f"Shaded detectors are used\n in this observation",
            color='cyan', transform=ax.transAxes, horizontalalignment='right')
    plt.text(0.98, 0.02, f"Cyan = Science attitude",
            color='cyan', transform=ax.transAxes, horizontalalignment='right')


    # re-establish the axes limits
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    return ax

##--- Functions for plotting the visit field of REGARD

def plot_circle_coords(x1, x2, from_frame, to_frame, ax, **plot_kwargs):
    """Plot a circle (great or otherwise) with astronomical coordinate transforms"""
    # make a circle of coords
    from_plane = coords.SkyCoord(x1, x2, unit='deg', frame=from_frame)
    return plot_coords_in_frame(from_plane, ax, to_frame, **plot_kwargs)

def plot_coords_in_frame(skycoords, ax, to_frame, **plot_kwargs):
    """ Plot in ICRS or Barycentric coordinates, with appropriate transformations
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
    ax.plot(-plot_x, plot_y, **plot_kwargs)



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

    #plt.text(0.99, 0.01, f"JWST field of regard\non {datetime.to_value('iso',subfmt='date')}\n[Ecliptic coords]",
    #         transform=ax.transAxes, horizontalalignment='right')

    plt.title(f"JWST field of regard\non {datetime.to_value('iso',subfmt='date')}\n[Ecliptic coords]",)
    
    # Plot all the markers
    plot_celestial_markers(visit, ax, 'geocentricmeanecliptic', show_sun=show_sun, datetime=datetime)

    return ax




def show_field_of_regard_ra_dec(visit, datetime=None, subplotspec=None, labelaxes=False):
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
    plot_celestial_markers(visit, ax, 'gcrs', datetime=datetime)

    return ax


def plot_celestial_markers(visit, ax, frame='gcrs', datetime=None, show_sun=True):
    """Main routine to plot celestial markers for sun, target, and various circles
    """

    if datetime is None:
        # What date/time are we generating this plot for?
        datetime = visit.time_early # astropy.time.Time.now()

    # Determine where the sun will be then.
    # NOTE - this is the apparent position as seen from Earth, and may be off by up to ~0.3 deg from the
    # perspective of JWST in its large orbit around L2. This is not a high accuracy calculation.
    sun = coords.get_sun(datetime).transform_to(frame)

    # Draw the sun, and anti-sun point
    antisun = sun.directional_offset_by(0, 180*u.deg)
    if show_sun: # it's less useful to show this in the lambert ecliptic proj
        plot_coords_in_frame(sun, ax, frame,  marker='o', markersize=20, color='orange', markeredgecolor='orange', zorder=10)
    plot_coords_in_frame(antisun, ax, frame,  marker='+', markersize=10, color='black', markeredgewidth=3,)

    # Draw galactic and ecliptic planes
    # JWST's CVZs are defined by the Barycentric version of the ecliptic plane.
    plot_circle_coords(np.arange(361), np.zeros(361), 'galactic', frame, ax,
                       ls='none',marker='.', markersize=2, color='maroon', alpha=0.3)
    plot_circle_coords(np.arange(361)-180, np.zeros(361), 'barycentricmeanecliptic', frame, ax,
                       ls='none',marker='.', markersize=2, color='black', alpha=0.3)
    # Draw the CVZs
    for ecliptic_lat in [85, -85]:
        plot_circle_coords(np.arange(361) - 180, np.zeros(361)+ecliptic_lat, 'barycentricmeanecliptic', frame, ax,
                           ls='none', marker='.', markersize=0.5, color='blue', alpha=0.3)


    # Shade the allowed field of regard.
    # Given hassles with plotting filled shapes in map projections, a relatively straightforward way to do this
    # ends up being to rasterize a grid over the sky, and scatter plot those points densely
    n = 180
    lat, lon = np.meshgrid(np.linspace(-np.pi, np.pi, 2 * n), np.linspace(-np.pi / 2, np.pi / 2, n))
    skymesh = astropy.coordinates.SkyCoord(lat, lon, frame=frame, unit='rad')
    #field_of_regard_mask = (seps > 85 * astropy.units.deg) & (seps < 135 * astropy.units.deg)

    pas = np.arange(360) * u.deg
    for offset in [85, 145]:
        for_edge = sun.directional_offset_by(pas, offset*u.deg)
        plot_coords_in_frame(for_edge, ax, frame,  marker='.',color='green', markersize=1, ls='none')
        
    for offset in np.linspace(86, 144, 50):
        for_circ = sun.directional_offset_by(pas, offset * u.deg)
        plot_coords_in_frame(for_circ, ax, frame, marker='o', color='#E5F2E5', markersize=1, ls='none', zorder=-30)
    #
    # if 0:
    #     with warnings.catch_warnings():
    #         # Temporarily ignore any benign warnings of math errors in the following call
    #         warnings.simplefilter("ignore")
    #         # Note, must flip sight on X coord given matplotlib convention vs. RA
    #         if 'ecliptic' in frame:
    #             plt.scatter(-lat[field_of_regard_mask], lon[field_of_regard_mask], color='#E5F2E5', zorder=-1000);
    #         else:
    #             plt.scatter(-skymesh[field_of_regard_mask].ra.wrap_at('180d').radian,
    #                         skymesh[field_of_regard_mask].dec.radian, color='#E5F2E5', zorder=-1000);

    # Show the target pointing!
    gs = coords.SkyCoord(visit.slew.GSRA, visit.slew.GSDEC, unit='deg', frame='icrs')
    plot_coords_in_frame(gs, ax, frame,  marker='*',color='red', markersize=20,)

    gal_center = astropy.coordinates.SkyCoord(0, 0, unit='deg', frame='galactic').transform_to('icrs')
    plot_coords_in_frame(gal_center, ax, frame, marker='o', markersize=5, color='maroon')

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
    fgs_aperture = visit.get_guider_aperture()
    fgs_aperture.set_attitude_matrix(attmat)
    #v1coords = astropy.coordinates.SkyCoord(*fgs_aperture.tel_to_sky(0,0), unit='deg', frame='icrs')
    v1ra, v1dec = fgs_aperture.tel_to_sky(0,0)
    v1ra = np.deg2rad(v1ra)
    v1dec = np.deg2rad(v1dec)

    #sun = astropy.coordinates.get_sun(times)

    # Load JWST ephemeris from jwst_gtvt
    A_eph = EPH.Ephemeris(
        os.path.join(os.path.dirname(os.path.abspath(jwst_gtvt.__file__)), "horizons_EM_jwst_wrt_sun_2020-2024.txt"),
        False, verbose=False)

    v3pa_nominal = np.zeros(nt, float)
    max_allowed_roll = np.zeros(nt, float)
    sun_angle = np.zeros(nt, float)


    for i in range(nt):
        # the following all in radians, for compatibility with jwst_gtvt internals
        atime = times[i]
        v3pa_nominal[i] = A_eph.normal_pa(atime.mjd, v1ra, v1dec)
        # Determine apparent sun position at that time.
        sun_ra, sun_dec = A_eph.sun_pos(atime.mjd)

        sun_angle[i] = jwst_gtvt.find_tgt_info.angular_sep(sun_ra, sun_dec, v1ra, v1dec)

        max_allowed_roll[i] = jwst_gtvt.find_tgt_info.allowed_max_vehicle_roll(sun_ra, sun_dec, v1ra, v1dec)

    # convert to degrees
    sun_angle = np.rad2deg((sun_angle))
    max_allowed_roll = np.rad2deg((max_allowed_roll))
    v3pa_nominal = np.rad2deg((v3pa_nominal))

    vehicle_pitch = sun_angle - 90

    V3PA = visit.slew.GSPA - 1.25
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
            ax.legend(loc='upper right')

    if np.any( np.abs(vehicle_roll) > max_allowed_roll):
        plt.text(0.98, 0.1, 'WARNING! Vehicle roll is outside of limits at some times!', color='red',
                 fontweight='bold', horizontalalignment='right',fontsize=15,
                 transform=plt.gcf().transFigure)


def multi_plot(visit, verbose=False, save=False):
    """ Main top-level function for visitviewer"""
    fig = plt.figure(figsize=(16, 9))

    # Set up a bunch of plot axes, complicatedly via nested grids
    # (this could be simplified but I'm adapting from other existing code...)
    gs_outer = plt.GridSpec(1, 2, width_ratios=[9,7])
    gs_outer.update(bottom=0.05, top=0.95, wspace=0.2, left=0.02, right=0.98)

    r_gridspec_kw = {'hspace': 0.2, 'wspace': 0.2, 'height_ratios': [0.3, 1, 1, 0.3], 'width_ratios': [1, 1.5]}
    gs_r = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=gs_outer[1],
                                            **r_gridspec_kw)
    r2_gridspec_kw = {'hspace': 0.05, 'wspace': 0.2, 'height_ratios': [0.3, .7, 0.3, 0.3, 0.4]}
    gs_r2 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs_outer[1],
                                             **r2_gridspec_kw)

    # Now make the plots, via the functions defined above
    plot_visit_fov(visit, subplotspec=gs_outer[0])
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
        plt.savefig(outname)
        print(f"File saved to {outname}")