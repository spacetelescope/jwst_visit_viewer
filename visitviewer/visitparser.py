import re
import os
import astropy
import pysiaf
import warnings

#######  Visit Quick Viewer
#
#  Code to make a simple, quick, human-friendly view of a visit file and its contents
#
#
######


# --- Part 1: Load and parse visit file information


iswhitespace = lambda x: re.fullmatch("\\s+", x) is not None


# Some lightweight classes for parsed information
# We use these to hold information in a convenient set of nested structures, and print it neatly

class Statement(object):
    """Any kind of statement in a visit file"""

    def __init__(self, cmdstring, verbose=False):

        cmdparts = cmdstring.split(' ,')
        self.name = cmdparts[0]
        self.args = cmdparts[1:]

        for param in self.args:
            if '=' in param:
                key, value = param.split('=')
                try:
                    value = float(value)
                except ValueError:
                    pass
                self.__dict__[key] = value

        if verbose:
            for part in cmdparts:
                print("   " + part)

    def __repr__(self):
        return "<Statement {}: {} >".format(self.name, self.args)


class SlewOrActStatement(Statement):
    """Base class for either a slew or activity statement, within some group in a visit file

    """

    def __init__(self, cmdstring, group=None, sequence=None, verbose=False):
        super().__init__(cmdstring, verbose=verbose)
        self.group = group
        self.sequence = sequence
        self.activity = int(self.args[0], 36)  # note, these are base 36 numbers; each digit counts 0 to 9 then A to Z
        self.scriptname = self.args[1]

        for param in self.args[2:]:
            if '=' in param:
                key, value = param.split('=')
                try:
                    value = float(value)
                except ValueError:
                    pass
            self.__dict__[key] = value

    @property
    def gsa(self):
        """Group, sequence, activity"""
        return "{:02d}{:1d}{:02d}".format(self.group, self.sequence, self.activity)


class SlewStatement(SlewOrActStatement):
    """Slew statement"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # OSS defaults to Fine Guide if GUIDEMODE is omitted
        # see https://innerspace.stsci.edu/display/OPGS/OSS+8.4+%28Forms+7.2%29+Pointing+OPGS+Rules#OSS8.4(Forms7.2)PointingOPGSRules-Slew
        if not hasattr(self, "GUIDEMODE"):
            self.GUIDEMODE = 'FINEGUIDE'

        # some slew statements are not actually slews!
        self._is_no_slew = self.scriptname == 'SCNOSLEWMAIN'

    def __repr__(self):
        if self._is_no_slew:
            return f"Slew {self.gsa}: SCNOSLEWMAIN"
        else:
            return "Slew {}: for {} on GS at ({}, {}) with GSPA={}".format(self.gsa, self.GUIDEMODE, self.GSRA, self.GSDEC,
                                                                        self.GSPA)


class ActivityStatement(SlewOrActStatement):
    """Activity statement, typically a call to an OSS script"""

    def __init__(self, cmdstring, *args, **kwargs):
        super().__init__(cmdstring, *args, **kwargs)

    def __repr__(self):
        return "Activity {}:  {}".format(self.gsa, self.describe())

    def describe(self):
        if self.scriptname == 'NRCWFSCMAIN':
            description = """{s.scriptname}  {s.CONFIG} WFCGROUP={s.WFCGROUP}
        Readout={s.NGROUPS:.0f} groups, {s.NINTS:.0f} ints
        SW={s.FILTSHORTA}, LW={s.FILTLONGA}"""
        elif self.scriptname == 'NRCMAIN' and hasattr(self, 'CONFIG'):
            module = 'B' if 'B' in self.CONFIG else 'A'
            description = """{s.scriptname}  {s.CONFIG}
        Readout={s.NGROUPS:.0f} groups, {s.NINTS:.0f} ints
        SW={s.FILTSHORT""" + module + "}, LW={s.FILTLONG" + module+"}"
        elif self.scriptname == 'NRCMAIN' and hasattr(self, 'CORON'):
            module = 'A'
            description = """{s.scriptname}  coron={s.CORON}
        Readout={s.NGROUPS:.0f} groups, {s.NINTS:.0f} ints
        filter={s.CORFILT}"""
        elif self.scriptname == 'NRCWFCPMAIN':
            if hasattr(self, "WFCGROUP"):
                description = f"{self.scriptname}  mirror move WFCGROUP={int(self.WFCGROUP)}"
            else:
                # DHS exposures
                mod = self.CONFIG[3]  # a or B
                description = "{s.scriptname}  NRC" + mod + ", {s.FILTSHORT" + mod + "}+{s.PUPILSHORT" + mod + \
                              "}, Readout={s.NGROUPS:.0f} groups, {s.NINTS:.0f} ints"
        elif self.scriptname == 'SCSAMMAIN':
            # SAMs can be set via various possible sets of parameters
            if hasattr(self, 'DELTAX'):
                description = """SCSAMMAIN  dx={s.DELTAX}, dy={s.DELTAY}, dpa={s.DELTAPA}"""
            elif hasattr(self, 'FGS1DELTAX'):
                description = """SCSAMMAIN  FGS1_deltax={s.FGS1DELTAX}, FGS1_dy={s.FGS1DELTAY}, dpa={s.FGS1DELTAPA}"""
            else:
                # to be improved if we encounter any other formats
                description = """SCSAMMAIN"""


        elif self.scriptname == 'NRCSUBMAIN':
            description = """NRCSUBMAIN   subarray={s.SUBARRAY}"""
        else:
            description = """{s.scriptname}"""
        return description.format(s=self)


class GuideStatement(SlewOrActStatement):
    """Guide statement"""

    def describe(self):
        if self.args[1] == 'FGSVERMAIN':
            return "Verification"
        else:
            detnum = self.DETECTOR[-1]
            return f"""FGS{detnum}\n        {", ".join(self.args[2:])}"""

    def __repr__(self):
        return "Guide {}:  {}".format(self.gsa, self.describe())


def parse_visit_file(lines):
    # Read the APT template name from the first line as a comment
    if lines[0].startswith("# "):
        apt_template = lines[0][2:].strip()
    else:
        apt_template = "UNKNOWN"

    # for parallels, second line will give a comment
    if lines[1].startswith("# ") and not lines[1].startswith("# OSS"):
        apt_template += " with parallel "+lines[1][2:].strip()


    # Simple parsing that ignores commands and newlines, but respects the fact that
    # OSS parameters are separated by the exact string " ," with the comma necessarily after whitespace.
    nocomments = [line.strip() for line in lines if not (line.startswith("#") or iswhitespace(line))]
    for i in range(len(nocomments)):
        if nocomments[i].startswith(','):
            nocomments[i] = ' ' + nocomments[i]

    merged = ''.join(nocomments)
    commands = merged.split(';')

    # Now iterate through the statements
    statements = []
    ct_group = 0
    ct_seq = 0

    for cmd in commands:
        if cmd == '':
            continue
        parsedcmd = Statement(cmd)
        if parsedcmd.name == 'GROUP':
            ct_group = int(parsedcmd.args[0].split()[0])
            continue  # No need to save this as separate statement?
        elif parsedcmd.name == 'SEQ':
            ct_seq = int(parsedcmd.args[0].split()[0])
            continue  # No need to save this as separate statement?
        elif parsedcmd.name == 'SLEW':
            parsedcmd = SlewStatement(cmd, group=ct_group, sequence=ct_seq)
        elif parsedcmd.name == 'ACT':
            if parsedcmd.args[1] == 'FGSMAIN' or parsedcmd.args[1] == 'FGSVERMAIN' or  parsedcmd.args[1] == 'FGSMTMAIN' :
                parsedcmd = GuideStatement(cmd, group=ct_group, sequence=ct_seq)
            else:
                parsedcmd = ActivityStatement(cmd, group=ct_group, sequence=ct_seq)
        statements.append(parsedcmd)
    return statements, apt_template


class VisitFileContents(object):
    """ Wrapper object for visit file contents

    Parse a visit file into a set of nested objects we can work with more conveniently.

    """

    def __init__(self, filename, no_gspa_yoffset=False, verbose=True, **kwargs):
        if verbose: print(f"Parsing {filename}")

        # Find the basics
        self.filename = filename
        with open(filename) as file:
            self.text = [line.rstrip() for line in file.readlines()]

        self.statements, self.template = parse_visit_file(self.text)
        self.visitid = os.path.splitext(os.path.basename(self.filename))[0]

        # Find the visit setup precursors
        self.dither = self._find_statement('DITHER')
        self.aux = self._find_statement('AUX')
        self.momentum = self._find_statement('MOMENTUM')

        # Find visit timing info and parse into astropy Times
        self.visit_timing = self._find_statement('VISIT')
        self.time_early = astropy.time.Time.strptime(self.visit_timing.EARLY, '%Y-%j/%H:%M:%S')
        self.time_late = astropy.time.Time.strptime(self.visit_timing.LATE, '%Y-%j/%H:%M:%S')
        # CUTOFF needs extra handling for :CONVST appendage
        # self.time_cutoff = astropy.time.Time.strptime(self.visit_timing.CUTOFF, '%Y-%j/%H:%M:%S')

        # Now find the actual visit activities
        # there is usually only one slew but sometimes there may be multiple
        self.slews = self._find_statements('SLEW')
        self.slew = self.slews[0]
        self.activities = self._find_statements('ACT')  # all the activities


        self.guide_activities = [a for a in self.activities if (isinstance(a, GuideStatement) and
                                 (a.scriptname != 'FGSVERMAIN'))]
        self.si_activities = [a for a in self.activities if (not isinstance(a, GuideStatement) and
                                                             (a.scriptname != 'GENWAITMAIN'))]

        # Record whether or not to interpret the GSPA parameter as having an angular offset of FGS V3IdlYangle, or not
        # PPS 14.14.1 and earlier do not apply this offset, but later versions more correctly do so.
        # This affects the inferred attitude matrix.
        self._no_gspa_yoffset = no_gspa_yoffset

        self.verbose = verbose

    def _find_statements(self, name):
        """ find all matching statements, return as list"""
        return [s for s in self.statements if s.name == name]

    def _find_statement(self, name):
        """ find all matching statements; return as scalar if there's just one, else return list"""
        matches = self._find_statements(name)
        if len(matches) == 1:
            return matches[0]
        return matches

    def _activity_apertures_used(self, activity):
        """ What are the main apertures used in this activity, which we should highlight on the plot?

        This function checks just a single OSS activity; see apertures_used() for the full version
        across the entire visit.

        Return list of SIAF aperture names
        """
        apernames = []
        if self.verbose: print(activity.scriptname)
        if activity.scriptname.startswith("NRC"):
            if activity.scriptname != 'NRCSUBMAIN':
                # NIRCam full-frame readouts, one or both modules
                if activity.scriptname == 'NRCIPRMAIN':
                    config = 'NRCALL'  # special case for the IPR template which does not have a 'CONFIG' parameter
                elif activity.scriptname == 'NRCPILMAIN' :
                    config = None  # This moves the PIL but does not take images
                elif activity.scriptname == 'NRCTAMAIN' :
                    config = None  # This performs TA, but we can ignore that and read the apertures from some later NRCMAIN call
                elif activity.scriptname == 'NRCWFCPMAIN' :
                    try:
                        config = activity.CONFIG
                    except AttributeError:
                        config = None  # special case for Coarse Phasing, which sometimes moves mirrors instead of taking images. 
                                       # thos calls do not have a 'CONFIG' parameter
                elif activity.scriptname == 'NRCMAIN' and hasattr(activity, 'CORON'):
                    config = 'CORON'
                    coron_aper_map = {'MASKA210R': 'NRCA2_MASK210R',
                                      'MASKA335R': 'NRCA5_MASK335R',
                                      'MASKA430R': 'NRCA5_MASK430R',
                                      'MASKASWB': 'NRCA4_MASKSWB',
                                      'MASKALWB': 'NRCA5_MASKLWB',
                            }
                    apernames.append(coron_aper_map[activity.CORON])
                else:
                    config = activity.CONFIG

                if self.verbose: print(f"Activity: {activity.scriptname}\tConfig: {config}")
                if config == 'NRCAALL' or config == 'NRCALL' or config == 'NRCASHORT':
                    [apernames.append(f'NRCA{n}_FULL') for n in (1, 2, 3, 4)]
                if config == 'NRCBALL' or config == 'NRCALL' or config == 'NRCBSHORT':
                    [apernames.append(f'NRCB{n}_FULL') for n in (1, 2, 3, 4)]
            else:
                # NRCSUBMAIN configures NIRCam subarray readouts
                if activity.SUBARRAY=='SUB96DHSPILA':
                    apernames.append('NRCA3_DHSPIL_SUB96')
                elif activity.SUBARRAY == 'SUB96DHSPILB':
                    apernames.append('NRCB4_DHSPIL_SUB96')
                # Add other NIRCam subarrays of interest here...

        elif activity.scriptname == 'SCWFCMAIN':
            # WF control 'move mirrors' script; no images taken
            pass
        elif activity.scriptname == 'SCSAMMAIN':
            # spacecraft small angle maneuver script; no images taken
            pass
        elif activity.scriptname == 'FGSIMAGEMAIN':
            # fgs_siaf = pysiaf.Siaf('FGS')
            # DETECTOR is 'GUIDER1' or 'GUIDER2'; convert to OSS aperture name
            apernames.append(f'FGS{activity.DETECTOR[6]}_FULL')
        # SI TEAMS PLEASE ADD IMPROVED ELIF STATEMENTS HERE FOR HOW TO INFER APERTURES FOR YOUR VARIOUS TEMPLATES
        # IF YOU WOULD LIKE MORE PRECISE PLOTS THAN THESE DEFAULTS
        elif activity.scriptname.startswith('NIS'):
            apernames.append("NIS_CEN")
        elif activity.scriptname.startswith('MIR'):
            apernames.append("MIRIM_FULL")
        elif activity.scriptname.startswith('NRS'):
            for n in range(4):
                apernames.append('NRS_FULL_MSA' + str(n + 1))
        else:
            warnings.warn(f"Inference of main apertures not implemented for that script, {activity.scriptname}...")

        return apernames

    def apertures_used(self):
        """ which SIAF apertures are the main ones used in this for data collection?
        This is a loose definition, from the perspective of the WF team currently.

        This is used in visitviewer to determine which apertures to shade in the plot.

        """

        main_apertures = []
        for activity in self.si_activities:
            main_apertures += self._activity_apertures_used(activity)

        # If we're guiding, that aperture counts too
        if self.slew.GUIDEMODE != 'COARSE':
            for activity in self.guide_activities:
                main_apertures.append(f"FGS{activity.DETECTOR[-1]}_FULL_OSS")

        # eliminate duplicates
        main_apertures = list(set(main_apertures))

        return main_apertures

    def uses_guiding(self):
        return len(self.guide_activities) > 0

    def uses_sams(self):
        return 'SCSAMMAIN' in [a.scriptname for a in self.si_activities]

    def short_summary(self):
        """ Succinct summary of what's in this visit. Should be no more than 1-2 lines

        Used in visitviewer for display at upper right.
        """
        text = self.template
        if self.template.startswith("WFSC"):
            act = self.si_activities[0]
            if self.template.startswith("WFSC NIRCam Coarse Phasing") and act.scriptname=='NRCSUBMAIN':  # skip over the subarray calls at the start of a coarse phasing WFSC visit
                act = self.si_activities[3]
            #text += f"\n    using {nrc_module}"
            mod = act.CONFIG[3]  # a or B
            description = "NRC " + mod + ", {act.FILTSHORT" + mod + "}+{act.PUPILSHORT" + mod + "}" +\
                          ", readout {act.PATTERN}, NGROUPS={act.NGROUPS:.0f}, NINTS={act.NINTS:.0f}"
            try:
                text += "\n  using "+ description.format(act=act)
            except AttributeError:
                pass # if we can't format the above just continue benignly

        return text

    def summarize(self):
        """Longer summary of visit file information

        Note, this was written before VPR, and is much simpler / less full-featured than that.
        """

        print("==Summary for file {}==".format(self.filename))
        print("From APT template: {}".format(self.template))
        # Check slew statements
        slews = [s for s in self.statements if (isinstance(s, Slew) or isinstance(s, Guide))]
        # Check activity statements
        acts = [s for s in self.statements if isinstance(s, Activity)]
        is_wfsc_visit = any(['WFSC' in a.scriptname for a in acts])

        print("\nContains {} slew or guide statement(s):".format(len(slews)))
        for s in slews:
            print("   " + repr(s))
        print("\nContains {} activity statement(s):".format(len(acts)))
        act_by_num = dict()
        for s in acts:
            print("   " + repr(s))
            act_by_num[s.gsa] = s
        if is_wfsc_visit:
            aux = [s for s in self.statements if s.name == 'AUX']
            if len(aux) == 0:
                raise RuntimeError("WFSC VISIT BUT NO AUX STATEMENT FOUND!")
            # Check for presence of AUX statement

    def get_attitude_matrix(self, step='slew', gscandidate=1, fgs_delta_from_sam=None):
        """Return attitude matrix for 'id' or 'science' attitudes in a visit

        For reference, see https://innerspace.stsci.edu/display/OPGS/OSS+8.4+%28Forms+7.2%29+Pointing+OPGS+Rules#OSS8.4(Forms7.2)PointingOPGSRules-Slew

        For FINEGUIDE or TRACK there is always a call to SCSLEWMAIN followed by a call to FGSMAIN, with similar parameters
        provided to both. For convenience we pull these out from the SLEW statement, at least for now.

        TODO: improve/enhance for more complicated visits with multiple guide stars, multiple slews, etc?

        Parameters
        -----------
        step : basestring
            Whether to read from slew statement (e.g. slew to ID attitude), or from science attitude
        fgs_delta_from_sam : None, or 2-tuple of floats
            Offset to apply in the FGS frame, as from an offset Small Angle Maneuver via SCSAMMAIN.
            Interpreted as (FGS1DELTAX, FGS1DELTAY) offset provided in FGS Ideal coordinates in units of arcseconds
        gscandidate: int
            Which GS candidate to draw the informayion from. 1-based index
        """
        visit = self

        indx = gscandidate-1

        # Read pointing information from SLEW or FGSMAIN statements in the visit
        if step == 'slew':
            ra = visit.slews[indx].GSRA
            dec = visit.slews[indx].GSDEC
            pa = visit.slews[indx].GSPA
            x_idl = visit.slews[indx].GSX
            y_idl = visit.slews[indx].GSY
        elif step == 'id':
            if not self.uses_guiding():
                raise RuntimeError("This visit does not use guiding, therefore there is no Guide ID attitude.")
            # is this always the same as the slews[indx]?

            # print(step, indx)
            try:
                ra = visit.guide_activities[indx].GSRA
                dec = visit.guide_activities[indx].GSDEC
            except AttributeError:
                # subsequent GS calls (not the first in a visit file) may not have GSRA, GSDEC in the FGSMAIN call,
                # so look at the previous slew statement?
                ra = visit.slews[indx].GSRA
                dec = visit.slews[indx].GSDEC

            pa = visit.guide_activities[indx].GSROLLID
            x_idl = visit.guide_activities[indx].GSXID
            y_idl = visit.guide_activities[indx].GSYID
            #print(ra,dec, pa, x_idl, y_idl, visit.guide_activities[indx].gsa)
        elif step == 'sci':
            if self.uses_guiding():
                # Guided visit, therefore read science attitude from the FGSMAIN call
                try:
                    ra = visit.guide_activities[indx].GSRA
                    dec = visit.guide_activities[indx].GSDEC
                except AttributeError:
                    # subsequent GS calls (not the first in a visit file) may not have GSRA, GSDEC in the FGSMAIN call,
                    # so look at the previous slew statement?
                    ra = visit.slews[indx].GSRA
                    dec = visit.slews[indx].GSDEC
                try:
                    pa = visit.guide_activities[indx].GSROLLSCI
                    x_idl = visit.guide_activities[indx].GSXSCI
                    y_idl = visit.guide_activities[indx].GSYSCI
                except AttributeError:
                    # Fall back to these if there is no science attitude provided?
                    pa = visit.guide_activities[indx].GSROLLID
                    x_idl = visit.guide_activities[indx].GSXID
                    y_idl = visit.guide_activities[indx].GSYID
            else:
                # Unguided visit, therefore infer science attitude must be same as slew attitude.
                ra = visit.slews[indx].GSRA
                dec = visit.slews[indx].GSDEC
                pa = visit.slews[indx].GSPA
                x_idl = visit.slews[indx].GSX
                y_idl = visit.slews[indx].GSY
        else:
            raise ValueError("step must be one of {slew, id, sci}.")

        fgs_aperture = self.get_guider_aperture()

        if fgs_delta_from_sam is not None:
            # Offset the guide star for a SAM. We interpret this as changing which location in the FGS frame
            # should be pointed toward the guide star or other target ra,dec
            if self.verbose:
                print(f"Applying SAM offset {fgs_delta_from_sam} in arcsec, FGS Ideal frame")
            x_idl = x_idl + fgs_delta_from_sam[0]
            y_idl = y_idl + fgs_delta_from_sam[1]


        # Convert from FGS1/2 Ideal coords in arcseconds to telescope FOV coordinates in arcseconds
        xtel, ytel = fgs_aperture.idl_to_tel(x_idl, y_idl)

        # The spacecraft ACS is expecting FGS Yics PA, which may or may not be what's in the visit file.
        # Here we ourselves want to get the V3PA, because that's what the attitude_matrix transform function wants.
        if self._no_gspa_yoffset:
            # The provided GSPA should be interpreted directly as the V3 PA at that guide star
            v3pa_at_guidestar = pa
        else:
            # The provided GSPA should be interpreted as the PA of the FGS Ideal coordinate system Y axis,
            # which is rotated relative to the V3 axis by some amount
            # We learned in flight this is ALWAYS interpreted as in the FGS *1* frame, regardless of which FGS is being used


            #fgs_Yics_offset = fgs_aperture.V3IdlYAngle  # Degrees, rotation offset between FGS1 Yics and V3PA, from SIAF
            #fgs_Yics_offset = -1.25
            from .visitplotter import SIAFS
            fgs_Yics_offset = SIAFS['FGS'].apertures['FGS1_FULL_OSS'].V3IdlYAngle
            v3pa_at_guidestar =  pa - fgs_Yics_offset

        # Compute attitude matrix
        attmat = pysiaf.utils.rotations.attitude_matrix(xtel, ytel,
                                                        ra, dec, v3pa_at_guidestar)
        return attmat

    def get_guider_aperture(self, return_name=False):
        """Return SIAF aperture for guider 1 or 2, as relevant
        parameters
        ----------
        return_name : bool
            default is to return the aperture object ; set this to
            True to return just the string aperture name
        """
        fgs_detector = 1 if self.slew.DETECTOR == 'GUIDER1' else 2
        fgs_aperture_name = f"FGS{fgs_detector}_FULL_OSS"
        if return_name:
            return fgs_aperture_name
        else:
            from .visitplotter import SIAFS
            fgs_aperture = SIAFS['FGS'].apertures[fgs_aperture_name]
            return fgs_aperture

    def get_nrcfilecount(self):
        """Return count of NIRCam files expected to be produced by the visit

        Contributed by Karl Misselt, NIRCam team, U of Arizona
        ----------
        """
        nrcFileCount = 0
        # might need others...
        # Don't really know what I'm doing here in detail, but should work for things like NRCMAIN and NRCSUBMAIN+NRCMAIN
        # NRCSUBMAIN is a subarray configuration call and MUST (?) be followed by an expose CONFIG
        # Don't know how to count NRCIPRMAIN, NRCCORON, NRCTSIMAGE need to find visit examples to see
        # what other tags to pull out.
        # AllowedScripts are _sciptnames; called in SUBMAIN+MAIN calls.
        # Some NRCWFCPMAIN calls don't seem to produce an obvious exposure in visit files...
        #    (MP replies: that's correct, those calls move mirrors rather than taking exposures)
        AllowedScripts = ["NRCMAIN", "NRCSUBMAIN", "NRCTAMAIN", "NRCFAMMAIN", "NRCWFCPMAIN", "NRCIPRMAIN", "NRCCORON",
                          "NRCTSIMAGE"]
        # AllowedConfigs define the commanded detector combinations. Depending on ExpType/TargetType, this will
        # define the number of SCA. DefaultNSCA is what we would naively expect
        AllowedConfigs = ["NRCA1", "NRCA2", "NRCA3", "NRCA4", "NRCALONG", "NRCASHORT", "NRCAALL",
                          "NRCB1", "NRCB2", "NRCB3", "NRCB4", "NRCBLONG", "NRCBSHORT", "NRCBALL",
                          "NRCALL"]
        DefaultNSca = [1, 1, 1, 1, 1, 4, 5, 1, 1, 1, 1, 1, 4, 5, 10]
        Config_To_NSca = dict(zip(AllowedConfigs, DefaultNSca))

        # May not need any of these defined?
        # AllowedExpTypes are "EXPTYPE=XXXX" in SUBMAIN calls
        # Why "FULL" lives in a SUBMAIN call, I've yet to grok.
        AllowedExpTypes = ["TA", "POINTSOURCE", "WF", "FULL", "LOS"]
        # AllowedTargetTypes are "TARGTYPE=XXXX" in NRCMAIN/NRCWFCPMAIN calls. And NRCSUBMAIN calls. Sheesh.
        # I *think* TIMESERIES is the only one I need worry about - nsca will always be 2 there, 1 SW, 1 LW
        AllowedTargetTypes = ["TIMESERIES", "EXTERNAL", "DARK"]

        for i in range(len(self.si_activities)):
            _scriptname = self.si_activities[i].scriptname
            if _scriptname in AllowedScripts:
                # For e.g. NRCSUBMAIN, there's no exposure taken, it's just configuring the ASIC.
                # There should be no other information we need from SUBMAIN since the configuration will be
                # repeated in the next activity. So for SUBMAIN we just skip to the next activity. Others?
                # TODO: Sort PIL imaging.
                #       Sort TIMESERIES
                #       Sort POINTSOURCE
                #       For POINTSOURCE, definition is in associated NRCSUBMAIN call; extract and carry
                #           over to next activity.
                #           Is POINTSOURCE (SUBMAIN) and TIMESERIES (MAIN) 1 or 2 detectors?
                #       How does GRISM get called? Always +3.  Definition is in NRCSUBMAIN call; extract
                #           and carry over to next activity

                # 'Null' _exposureType so that we only set it if the previous acticity was NRCSUBMAIN
                _exposureType = None
                if _scriptname == "NRCSUBMAIN":
                    # Save EXPTYPE here; continue and next statement should be NRCMAIN (Check/Verify algorithmically?)
                    _exposureType = self.si_activities[i].EXPTYPE
                    continue
                if _scriptname == "NRCTAMAIN":
                    # always a single SCA
                    nrcFileCount += 1
                elif _scriptname == "NRCCORON":
                    # always 1 scas
                    nrcFileCount += 1
                elif _scriptname == "NRCMAIN":
                    # Workhorse - Will need to check TIMESERIES, GRISM, POINTSOURCE, etc.
                    nrcFileCount += Config_To_NSca[self.si_activities[i].CONFIG]
                elif _scriptname == "NRCWFCPMAIN":
                    # Seems two types, one that generates exposures, and one that... doesn't; but does
                    # the mirror move. From jse script for expose, looks like 2 exposures per 'call'
                    # Differentiate between the 2 types of WFCP calls based on existence of WFCGROUP
                    # in the activity statement.
                    try:
                        self.si_activities[i].WFCGROUP
                        continue
                    except:
                        # WFCP generates 2 files per activity
                        _factor = 2
                        nrcFileCount += _factor * Config_To_NSca[self.si_activities[i].CONFIG]

        return nrcFileCount
