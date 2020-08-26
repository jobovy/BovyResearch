import matplotlib.pyplot as plt
import numpy
from astropy import units
from galpy.orbit import Orbit
from galpy.potential import ChandrasekharDynamicalFrictionForce
from galpy.potential import DehnenSmoothWrapperPotential
from galpy.potential import MWPotential2014 as mw
from galpy.potential import MovingObjectPotential
from galpy.util.bovy_coords import rect_to_cyl, rect_to_cyl_vec


def plot_satgalorbit(sat_galaxy, name_sat, xaxis_dim, yaxis_dim):

    # Function plots orbit (integrated backwards) of given satellite galaxy of the Milky Way Galaxy
    #
    # Input:
    # sat_galaxy = orbit object 
    # x and y axis dimensions to plot in 
    # name_sat: string of the satellite galaxy's name
    #
    # Output:
    #
    # Returns the plot of given orbit object
    
    ts = numpy.linspace(0.,-10.,1000)*units.Gyr
    sat_galaxy.integrate(ts, mw)
    plot = sat_galaxy.plot(d1 = xaxis_dim, d2 = yaxis_dim)
    plt.title('Orbit of Satellite Galaxy: '+ name_sat)
    plt.show()
    plt.close()
    return plot


def plot_satgalorbit_cdf(sat_galaxy, name_sat, xaxis_dim, yaxis_dim, sat_mass, sat_size):

    # Function plots orbits of satellite galaxy taking dynamical friction effects into account
    #
    # Input:
    # sat_galaxy = orbit object 
    # x and y axis dimensions to plot in (strings)
    # mass and size of satellite to model dynamical friction effects (quantities with units attached)
    # name_sat: string of the satellite galaxy's name
    #
    # Output:
    #
    # Returns the plot of given orbit object
    
    ts = numpy.linspace(0.,-10.,1000)*units.Gyr
    cdf= ChandrasekharDynamicalFrictionForce(GMs=sat_mass, rhm = sat_size ,dens=mw)
    sat_galaxy.integrate(ts, mw+cdf)
    plot = sat_galaxy.plot(d1 = xaxis_dim, d2 = yaxis_dim)
    plt.title('Orbit of Satellite Galaxy: '+ name_sat + ' Including Dynamical Friction')
    plt.show()
    plt.close()
    return plot


def plot_sat_cluster(sat_galaxy, name_sat, sat_potential, xaxis_dim, yaxis_dim, sat_mass, sat_size, x_satgal, y_satgal,
                     z_satgal, vx_satgal,vy_satgal, vz_satgal, tform, tsteady):

    # Function plots orbits of satellite galaxy as well as a star cluster within the satellite galaxy - simulates accretion onto MW
    #
    # Input:
    # sat_galaxy: an orbit object for a given satellite galaxy of the MW
    # name_sat: string of the satellite galaxy's name
    # sat_potential: potential object modelling the satellite's potential
    # x and y axis dimensions to plot in (strings)
    # x_satgal, y_satgal, z_satgal: x,y,z positions of the star_cluster within the satellite galaxy's frame of reference
    # vx_satgal, vy_satgal, vz_satgal: x,y,z velocities of the star_cluster within the satellite galaxy's frame of reference
    # mass and size of satellite to model dynamical friction effects (quantities with units attached)
    # tform, tsteady: parameters of the potential, models tidal disruption of satellite galaxy (quantities with units attached)
    #
    # Output:
    #
    # end_pos_cluster: the coordinates for star cluster at end of integration (@ time tend: tform+5*units.Gyr)
    # end_pos_gal: the coordinates for satellite galaxy at tend 
    # dswp: wrapper potential object that modifies satelltie galaxy to begin disrupting at tform 
    # cdf: model for dynamical friction force
     
    t_back = 10.
    ts = numpy.linspace(0.,-t_back,1000)*units.Gyr
    cdf= ChandrasekharDynamicalFrictionForce(GMs=sat_mass, rhm = sat_size ,dens=mw)
    sat_galaxy.integrate(ts, mw+cdf)
    '''
    R_sat = sat_galaxy.R(-t_back*units.Gyr)  #cylindrical radius at time t
    vR_sat = sat_galaxy.vR(-t_back*units.Gyr) #radial velocity at time t 
    vT_sat = sat_galaxy.vT(-t_back*units.Gyr) #tangential velocity at time t 
    z_sat = sat_galaxy.z(-t_back*units.Gyr) #vertical height at time t
    vz_sat = sat_galaxy.vz(-t_back*units.Gyr) #vertical velocity at time t 
    phi_sat = sat_galaxy.phi(-t_back*units.Gyr) #azimuth at time t 
    '''
    # Rectangular coordinates and velocities
    coord = [sat_galaxy.x(-t_back*units.Gyr), sat_galaxy.y(-t_back*units.Gyr), sat_galaxy.z(-t_back*units.Gyr)]
    vcoord = [sat_galaxy.vx(-t_back*units.Gyr),sat_galaxy.vy(-t_back*units.Gyr),sat_galaxy.vz(-t_back*units.Gyr)]
    
    t_fwrd = 15
    ts_f= numpy.linspace(-t_back, -t_back+t_fwrd,1000)*units.Gyr
    sat_galaxy = sat_galaxy(-t_back*units.Gyr)
    sat_galaxy.integrate(ts_f, mw + cdf)
    sat_galaxy.plot(d1 =  xaxis_dim, d2= yaxis_dim,linestyle = ':', color = 'black', label = 'satellite') #plots orbit of the satellite galaxy in MW frame of reference
    
    #!!sat_pot = HernquistPotential(amp = 2*sat_mass, a = sat_size, ro = 8., vo=220.)
    sat_movingpot = MovingObjectPotential(sat_galaxy, sat_potential)

    # Transform from satellite galaxy's frame of reference to Milky Way Galaxy's frame of reference (using Cartesian coordinates)
    # Rectangular coordinates of the star cluster in galactocentric frame
    x_gal = coord[0] + x_satgal
    y_gal = coord[1] + y_satgal
    z_gal = coord[2] + z_satgal
    # Velocity of the star cluster in galactocentric frame
    vx_gal = vcoord[0] + vx_satgal
    vy_gal = vcoord[1] + vy_satgal
    vz_gal = vcoord[2] + vz_satgal
    # Transform to cylindrical coordinate system: R, phi, z 
    R, phi, z = rect_to_cyl(x_gal, y_gal, z_gal)
    vR, vT, vz = rect_to_cyl_vec(vx_gal, vy_gal, vz_gal,x_gal, y_gal, z_gal, cyl = False)

    star_cluster = Orbit(vxvv = [R,vR,vT,z,vz,phi],ro = 8., vo=220.)
    star_cluster.integrate(ts_f, mw + sat_movingpot)
    star_cluster.plot(d1 =  xaxis_dim, d2= yaxis_dim, linestyle = '-', overplot = True, color = 'blue', alpha=0.6, label = 'star cluster') #plots orbit of the star_cluster in MW frame of reference
    plt.title('Orbit of Star Cluster Within Satellite Galaxy: ' + name_sat + ' in Galactocentric Frame')
    plt.legend()
    plt.show()
    plt.close()
    

    # Implement wrapper potential to simulate tidal disruption of satellite galaxy
    # Plot orbit of the satellite galaxy and star cluster within sat galaxy in MW frame of reference:
    plt.figure(figsize=(12.,10.))
    tstart = tform - 5.*units.Gyr
    tend = tform + 5.*units.Gyr 
    time_int = numpy.linspace(tstart.to_value(units.Gyr), tend.to_value(units.Gyr), 1000)*units.Gyr
    
    if tstart < -t_back*units.Gyr:
            # re-integrate satellite galaxy from current time back to tstart 
            re_time = numpy.linspace(-t_back, tstart.to_value(units.Gyr), 1000)*units.Gyr
            sat_galaxy.integrate(re_time, mw+cdf)
            
            # initialize star cluster on orbit in satellite galaxy at time tstart:
            # Rectangular coordinates and velocities
            coord = [sat_galaxy.x(tstart), sat_galaxy.y(tstart), sat_galaxy.z(tstart)]
            vcoord = [sat_galaxy.vx(tstart),sat_galaxy.vy(tstart),sat_galaxy.vz(tstart)]

            # Transform from satellite galaxy's frame of reference to Milky Way Galaxy's frame of reference (using Cartesian coordinates)
            # Rectangular coordinates of the star cluster in galactocentric frame
            x_gal = coord[0] + x_satgal
            y_gal = coord[1] + y_satgal
            z_gal = coord[2] + z_satgal
            # Velocity of the star cluster in galactocentric frame
            vx_gal = vcoord[0] + vx_satgal
            vy_gal = vcoord[1] + vy_satgal
            vz_gal = vcoord[2] + vz_satgal
            # Transform to cylindrical coordinate system: R, phi, z 
            R, phi, z = rect_to_cyl(x_gal, y_gal, z_gal)
            vR, vT, vz = rect_to_cyl_vec(vx_gal, vy_gal, vz_gal,x_gal, y_gal, z_gal, cyl = False)
            
            # Re-initialize star cluster on orbit at time tstart
            star_cluster = Orbit(vxvv=[R,vR,vT,vz,z,phi], ro=8., vo=220.)
    else:
        # default: star cluster is initialized at -10Gyr in given satellite galaxy
        star_cluster = star_cluster(tstart)
    
    
    sat_galaxy = sat_galaxy(tstart) #make copy of sat_galaxy orbit at time tstart 
    sat_galaxy.integrate(time_int, mw+cdf) # integrate sat_galaxy forward for 10Gyrs
    sat_galaxy.plot(d1 =  xaxis_dim, d2= yaxis_dim,linestyle = ':', color = 'black', label = 'satellite galaxy')
    sat_movingpot = MovingObjectPotential(sat_galaxy, sat_potential)
    dswp = DehnenSmoothWrapperPotential(amp=1.0, pot = sat_movingpot, tform=tform, tsteady=tsteady, decay = True)
    star_cluster.integrate(time_int, mw+dswp) 
        # star cluster in combined potential: MW galaxy & moving potential of satellite galaxy 
    star_cluster.plot(d1 =  xaxis_dim, d2= yaxis_dim, linestyle = '-', overplot = True, color = 'blue', alpha = 0.6,\
                      label = 'star cluster') 
        #plots orbit of the star_cluster in MW frame of reference
    plt.legend()
    plt.title('Star Cluster Orbit Within: '+name_sat+' for Tform = ' + str(tform) + ' & Tsteady = ' + str(tsteady) + ' in Galactocentric Frame')
    plt.savefig('WrapperPotential-Decaying Mass.pdf')
    plt.show()
    plt.close()
    
    # Figure out where star cluster is at end of integration: at tend 
    end_pos_cluster = [star_cluster.R(tend),star_cluster.vR(tend),star_cluster.vT(tend),star_cluster.z(tend),star_cluster.vz(tend), star_cluster.phi(tend)]
            # [R,vT,vT,z,vz,phi]
    end_pos_gal = [sat_galaxy.R(tend),sat_galaxy.vR(tend),sat_galaxy.vT(tend),sat_galaxy.z(tend),sat_galaxy.vz(tend), sat_galaxy.phi(tend)]

    '''
    # Used for finding dswp when integrating satellite galaxy backward in previous version of code
    time_intb = numpy.linspace(tend.to_value(units.Gyr), tstart.to_value(units.Gyr), 1000)*units.Gyr
    star_cluster_b = Orbit(vxvv = end_pos_cluster, ro=8., vo =220.) #full 6 coordinates
    sat_galaxy_b = Orbit(vxvv=end_pos_gal, ro=8., vo =220.)
    sat_galaxy_b.integrate(time_intb, mw + cdf)
    sat_galaxy_b.plot(d1 =  xaxis_dim, d2= yaxis_dim,linestyle = ':', color = 'black', label = 'satellite galaxy')
    sat_movingpot_b = MovingObjectPotential(sat_galaxy_b, sat_potential)
    #new_tform = tform - end_t
    #dswp_back = DehnenSmoothWrapperPotential(amp=1.0, pot = sat_movingpot_b, tform=tform, tsteady=tsteady, decay = True)
    star_cluster_b.integrate(time_intb, mw + dswp) # star cluster is in combined potential of MW galaxy and the moving potential of satellite galaxy 
    star_cluster_b.plot(d1 =  xaxis_dim, d2= yaxis_dim, linestyle = '-', overplot = True, color = 'blue', alpha = 0.6,\
                      label = 'star cluster') # galactocentric radius as a function of time
    plt.legend()
    plt.title('Orbit of Star Cluster Within Satellite Galaxy for Tform = ' + str(tform) + ' & Tsteady = ' + str(tsteady) + ' (in Galactocentric Frame)')
    plt.show()
    plt.close()
    '''
    

    
    return end_pos_cluster,end_pos_gal, dswp, cdf

