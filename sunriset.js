/*
 * sunriset.js
 * ===========
 *
 * Computes Sun rise/set and twilight times for any date and coordinates.
 *
 * Usage
 * -----
 * day_length(year, month, day, lon, lat)
 * day_civil_twilight_length(year, month, day, lon, lat)
 * day_nautical_twilight_length(year, month, day, lon, lat)
 * day_astronomical_twilight_length(year, month, day, lon, lat)
 * sun_rise_set(year, month, day, lon, lat)
 * civil_twilight(year, month, day, lon, lat)
 * nautical_twilight(year, month, day, lon, lat)
 * astronomical_twilight(year, month, day, lon, lat)
 *
 *
 * Example
 * -------
 * times = sun_rise_set(2009, 10, 14, 0.26, 51.56);
 * alert("Sunrise: " + times.trise);
 * alert("Sunset: " + times.tset);
 *
 *
 * Notes
 * -----
 * Longitudes and latitudes are specified in floating point degrees.
 * Longitudes EAST of Greenwich are POSITIVE.
 * Longitudes WEST of Greenwich are NEGATIVE.
 * Latitudes NORTH of the equator are POSITIVE.
 * Latitudes SOUTH of the equator are NEGATIVE.
 *
 * Months are numbered from 1 to 12 (so beware when using Date.getMonth() - you
 * will need to add 1.
 *
 * Years are expected in full 4-digit glory (so use Date.getFullYear() not
 * Date.getYear()).
 *
 * The day_* functions return a floating point hour value (e.g. 2.5 = 2h 30m.)
 *
 * The other four functions return an object with three properties:
 * rc: -1 if the sun is never above the horizon, +1 if the sun is always above
 *     the horizon, otherwise 0 (the sun rises and sets as normal.)
 * trise: time of sunrise (or start of period) as floating point hour value.
 * tset: time of sunset (or end of period) as floating point hour value.
 *
 *
 * Copyright and licence
 * ---------------------
 * Copyright (c) 2009, Neil de Carteret
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * sunriset.js isbased on SUNRISET.C by Paul Schlyter.
 * Copyright notice from SUNRISET.C:
 *
 * Written as DAYLEN.C, 1989-08-16
 * Modified to SUNRISET.C, 1992-12-01
 * (c) Paul Schlyter, 1989, 1992
 * Released to the public domain by Paul Schlyter, December 1992
 */


/*
 * Return the number of days elapsed since 2000 Jan 0.0
 * (which is equal to 1999 Dec 31, 0h UT)
 */
function days_since_2000_Jan_0 (y, m, d) {
    return 367*y - ((7*(y+((m+9)/12)))/4) + ((275*m)/9) +d-730530;
}

/*
 * Some conversion factors between radians and degrees
 */
var PI = Math.PI;
var RADEG = (180.0 / PI);
var DEGRAD = (PI / 180.0);
var INV360 = (1.0/360.0);

/*
 * The trigonometric functions in degrees
 */
function sind(x) { return Math.sin(x * DEGRAD); }
function cosd(x) { return Math.cos(x * DEGRAD); }
function tand(x) { return Math.tan(x * DEGRAD); }
function atand(x) { return RADEG * Math.atan(x); }
function asind(x) { return RADEG * Math.asin(x); }
function acosd(x) { return RADEG * Math.acos(x); }
function atan2d(y, x) { return RADEG * Math.atan2(y,x); }


/*
 * Following are some wrappers around the "workhorse" function __daylen__.
 * They mainly fill in the desired values for the reference altitude
 * below the horizon, and also selects whether this altitude should
 * refer to the Sun's center or its upper limb.
 */


/*
 * Compute the length of the day, from sunrise to sunset.
 * Sunrise/set is considered to occur when the Sun's upper limb is
 * 35 arc minutes below the horizon (this accounts for the refraction
 * of the Earth's atmosphere).
 */
function day_length(year, month, day, lon, lat) {
    return __daylen__(year, month, day, lon, lat, -35.0/60.0, 1);
}

/*
 * Compute the length of the day, including civil twilight.
 * Civil twilight starts/ends when the Sun's center is 6 degrees below
 * the horizon.
 */
function day_civil_twilight_length(year, month, day, lon, lat) {
    return __daylen__(year, month, day, lon, lat, -6.0, 0);
}

/*
 * Computes the length of the day, incl. nautical twilight.
 * Nautical twilight starts/ends when the Sun's center is 12 degrees
 * below the horizon.
 */
function day_nautical_twilight_length(year,month,day,lon,lat) {
    return __daylen__(year, month, day, lon, lat, -12.0, 0);
}

/*
 * Compute the length of the day, incl. astronomical twilight.
 * Astronomical twilight starts/ends when the Sun's center is 18 degrees
 * below the horizon.
 */
function day_astronomical_twilight_length(year,month,day,lon,lat) {
    return __daylen__(year, month, day, lon, lat, -18.0, 0);
}

/*
 * Computes times for sunrise/sunset.
 * Sunrise/set is considered to occur when the Sun's upper limb is
 * 35 arc minutes below the horizon (this accounts for the refraction
 * of the Earth's atmosphere).
 */
function sun_rise_set(year, month, day, lon, lat) {
    return __sunriset__(year, month, day, lon, lat, -35.0/60.0, 1);
}

/*
 * Compute the start and end times of civil twilight.
 * Civil twilight starts/ends when the Sun's center is 6 degrees below
 * the horizon.
 */
function civil_twilight(year, month, day, lon, lat) {
    return __sunriset__(year, month, day, lon, lat, -6.0, 0);
}

/*
 * Compute the start and end times of nautical twilight.
 * Nautical twilight starts/ends when the Sun's center is 12 degrees
 * below the horizon.
 */
function nautical_twilight(year, month, day, lon, lat) {
    return __sunriset__(year, month, day, lon, lat, -12.0, 0);
}

/*
 * Computes the start and end times of astronomical twilight.
 * Astronomical twilight starts/ends when the Sun's center is 18 degrees
 * below the horizon.
 */
function astronomical_twilight(year, month, day, lon, lat) {
    return __sunriset__(year, month, day, lon, lat, -18.0, 0);
}


/*
 * Calculate sun rise and set time for a given date, palce, and sun altitude.
 *
 * This is the workhorse function for sun rise/set times.
 *
 * ### Arguments
 * `year,month,date`: calendar date, 1801-2099 only.
 * `lon, lat`: longitude and latitude. See Notes.
 *      The longitude value IS critical in this function!
 * `altit`: the altitude which the Sun should cross.
 *      Set to -35/60 degrees for normal rise/set,
 *      -6 degrees for civil twilight, -12 degrees for nautical twilight, and
 *      -18 degrees for astronomical twilight.
 * `upper_limb`: non-zero -> upper limb, zero -> center.
 *      Set to non-zero (e.g. 1) when computing rise/set times, and to zero when
 *      computing start/end of twilight.
 *
 * ### Returns
 * An object with three properties:
 * `rc`: 0 = sun rises/sets this day. trise and tset are times.
 *      +1 = sun above the specified "horizon" 24 hours. trise set to time when
 *           the sun is at south, minus 12 hours while tset is set to the south
 *           time plus 12 hours. "Day" length = 24 hours.
 *      -1 = sun is below the specified "horizon" 24 hours.
 *           "Day" length = 0 hours, *trise and *tset are both set to the time
 *           when the sun is at south.
 * `trise`: time of sunrise (or start of period) as floating point hour value.
 * `tset`: time of sunset (or end of period) as floating point hour value.
 */
function __sunriset__(year, month, day, lon, lat, altit, upper_limb) {
    var  d,     // Days since 2000 Jan 0.0 (negative before)
    sr,         // Solar distance, astronomical units
    sRA,        // Sun's Right Ascension
    sdec,       // Sun's declination
    sradius,    // Sun's apparent radius
    t,          // Diurnal arc
    tsouth,     // Time when Sun is at south
    sidtime;    // Local sidereal time
    var ret = {};
    var rc = 0; // Return code from function - usually 0

    // Compute d of 12h local mean solar time
    d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

    // Compute local sideral time of this moment
    sidtime = revolution( GMST0(d) + 180.0 + lon );

    // Compute Sun's RA + Decl at this moment
    var my_sun_ra_dec = sun_RA_dec(d);
    sRA = my_sun_ra_dec.RA;
    sdec = my_sun_ra_dec.dec;
    sr = my_sun_ra_dec.r;

    // Compute time when Sun is at south - in hours UT
    tsouth = 12.0 - rev180(sidtime - sRA)/15.0;

    // Compute the Sun's apparent radius, degrees
    sradius = 0.2666 / sr;

    // Do correction to upper limb, if necessary
    if ( upper_limb )
        altit -= sradius;

    // Compute the diurnal arc that the Sun traverses to reach
    // the specified altitide altit:
    var cost;
    cost = ( sind(altit) - sind(lat) * sind(sdec) ) /
          ( cosd(lat) * cosd(sdec) );
    if ( cost >= 1.0 )
          rc = -1, t = 0.0;       // Sun always below altit
    else if ( cost <= -1.0 )
          rc = +1, t = 12.0;      // Sun always above altit
    else
          t = acosd(cost)/15.0;   // The diurnal arc, hours

    // Store rise and set times - in hours UT
    ret.trise = tsouth - t;
    ret.tset  = tsouth + t;
    ret.rc = rc;
    return ret;
}


/*
 * Calculate length of day for a given date, place, and sun altitude.
 *
 * This is the workhorse function for day lengths.
 *
 * ### Arguments
 * `year,month,date`: calendar date, 1801-2099 only.
 * `lon, lat`: longitude and latitude. See Notes.
 *      The longitude value is not critical in this function. Set it to the
 *      correct longitude if you're picky, otherwise set to to, say, 0.0.
 * `altit`: the altitude which the Sun should cross.
 *      Set to -35/60 degrees for normal rise/set,
 *      -6 degrees for civil twilight, -12 degrees for nautical twilight, and
 *      -18 degrees for astronomical twilight.
 * `upper_limb`: non-zero -> upper limb, zero -> center.
 *      Set to non-zero (e.g. 1) when computing rise/set times, and to zero when
 *      computing start/end of twilight.
 *
 * ### Returns
 * The length of the day in floating-point hours (e.g. 12h 30m is 12.5)
 */
function __daylen__(year, month, day, lon, lat, altit, upper_limb) {
    var  d,     // Days since 2000 Jan 0.0 (negative before)
    obl_ecl,    // Obliquity (inclination) of Earth's axis
    sr,         // Solar distance, astronomical units
    slon,       // True solar longitude
    sin_sdecl,  // Sine of Sun's declination
    cos_sdecl,  // Cosine of Sun's declination
    sradius,    // Sun's apparent radius
    t;          // Diurnal arc

    // Compute d of 12h local mean solar time
    d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

    // Compute obliquity of ecliptic (inclination of Earth's axis)
    obl_ecl = 23.4393 - 3.563E-7 * d;

    // Compute Sun's position
    var mysunpos = sunpos(d);
    lon = mysunpos.lon;
    sr = mysunpos.r;

    // Compute sine and cosine of Sun's declination
    sin_sdecl = sind(obl_ecl) * sind(slon);
    cos_sdecl = sqrt( 1.0 - sin_sdecl * sin_sdecl );

    // Compute the Sun's apparent radius, degrees
    sradius = 0.2666 / sr;

    // Do correction to upper limb, if necessary
    if ( upper_limb )
        altit -= sradius;

    // Compute the diurnal arc that the Sun traverses to reach
    // the specified altitide altit: */
    var cost;
    cost = ( sind(altit) - sind(lat) * sin_sdecl ) /
          ( cosd(lat) * cos_sdecl );
    if ( cost >= 1.0 )
          t = 0.0;                      // Sun always below altit
    else if ( cost <= -1.0 )
          t = 24.0;                     // Sun always above altit
    else  t = (2.0/15.0) * acosd(cost); // The diurnal arc, hours

    return t;
}


/*
 * Computes the Sun's position at any instant
 *
 * Computes the Sun's ecliptic longitude and distance
 * at an instant given in d, number of days since
 * 2000 Jan 0.0.  The Sun's ecliptic latitude is not
 * computed, since it's always very near 0.
 *
 * ### Arguments
 * d: number of days since 2000 Jan 0.0.
 *
 * ### Returns
 * An object with two properties:
 * lon: Sun's ecliptic longitude
 * r: Solar distance
 */
function sunpos(d) {
    var M,         // Mean anomaly of the Sun
        w,         // Mean longitude of perihelion
                   // Note: Sun's mean longitude = M + w
        e,         // Eccentricity of Earth's orbit
        E,         // Eccentric anomaly
        x, y,      // x, y coordinates in orbit
        v;         // True anomaly
    var pos = {};

    // Compute mean elements
    M = revolution( 356.0470 + 0.9856002585 * d );
    w = 282.9404 + 4.70935E-5 * d;
    e = 0.016709 - 1.151E-9 * d;

    // Compute true longitude and radius vector
    E = M + e * RADEG * sind(M) * ( 1.0 + e * cosd(M) );
        x = cosd(E) - e;
    y = Math.sqrt( 1.0 - e*e ) * sind(E);
    pos.r = Math.sqrt( x*x + y*y ); // Solar distance
    v = atan2d( y, x );             // True anomaly
    pos.lon = v + w;                // True solar longitude
    if ( pos.lon >= 360.0 )
        pos.lon -= 360.0;           // Make it 0..360 degrees
    return pos;
}


/*
 * Compute sun's right ascension and declination.
 *
 * ### Arguments
 * d: number of days since 2000 Jan 0.0.
 *
 * ### Returns
 * An object with two properties:
 * RA: Sun's right ascension.
 * dec: Sun's declination
 *
 * ### Notes
 * This function was undocumented in SUNRISET.C - I've filled it in with what I
 * can guess :)
 */
function sun_RA_dec (d) {
      var lon, obl_ecl, x, y, z, r;
      var ret = {};

      // Compute Sun's ecliptical coordinates
      mysunpos = sunpos(d);
      lon = mysunpos.lon;
      r = ret.r = mysunpos.r;

      // Compute ecliptic rectangular coordinates (z=0)
      x = r * cosd(lon);
      y = r * sind(lon);

      // Compute obliquity of ecliptic (inclination of Earth's axis)
      obl_ecl = 23.4393 - 3.563E-7 * d;

      // Convert to equatorial rectangular coordinates - x is uchanged
      z = y * sind(obl_ecl);
      y = y * cosd(obl_ecl);

      // Convert to spherical coordinates
      ret.RA = atan2d( y, x );
      ret.dec = atan2d( z, Math.sqrt(x*x + y*y) );

      return ret;
}


/*
 * Reduce angle to within 0..360 degrees
 */
function revolution (x) {
    return( x - 360.0 * Math.floor( x * INV360 ) );
}


/*
 * Reduce angle to within +180..+180 degrees
 */
function rev180 (x) {
    return( x - 360.0 * Math.floor( x * INV360 + 0.5 ) );
}


/*
 * Compute GMST0, the Greenwhich Mean Sidereal Time
 * at 0h UT (i.e. the sidereal time at the Greenwhich meridian at
 * 0h UT).  GMST is then the sidereal time at Greenwich at any
 * time of the day.  I've generelized GMST0 as well, and define it
 * as:  GMST0 = GMST - UT  --  this allows GMST0 to be computed at
 * other times than 0h UT as well.  While this sounds somewhat
 * contradictory, it is very practical:  instead of computing
 * GMST like:
 *
 *  GMST = (GMST0) + UT * (366.2422/365.2422)
 *
 * where (GMST0) is the GMST last time UT was 0 hours, one simply
 * computes:
 *
 *  GMST = GMST0 + UT
 *
 * where GMST0 is the GMST "at 0h UT" but at the current moment!
 * Defined in this way, GMST0 will increase with about 4 min a
 * day.  It also happens that GMST0 (in degrees, 1 hr = 15 degr)
 * is equal to the Sun's mean longitude plus/minus 180 degrees!
 * (if we neglect aberration, which amounts to 20 seconds of arc
 * or 1.33 seconds of time)
 */
function GMST0(d) {
    var sidtim0;
    /* Sidtime at 0h UT = L (Sun's mean longitude) + 180.0 degr  */
    /* L = M + w, as defined in sunpos().  Since I'm too lazy to */
    /* add these numbers, I'll let the C compiler do it for me.  */
    /* Any decent C compiler will add the constants at compile   */
    /* time, imposing no runtime or code overhead.               */
    sidtim0 = revolution( ( 180.0 + 356.0470 + 282.9404 ) +
                      ( 0.9856002585 + 4.70935E-5 ) * d );
    return sidtim0;
}
