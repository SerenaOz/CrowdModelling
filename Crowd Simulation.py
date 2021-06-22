from cromosim import *
from cromosim.micro import *
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from pylab import rcParams
rcParams['figure.figsize'] = 24.655, 3

def domainInit():
    width = 100.0
    height = 7.5

    # Create domain object
    dom = Domain(name = 'corridor',  pixel_size = 0.1, xmin=0, ymin=0,
                width=int(width*10), height=int(height*10))

    # Define wall and exit colours
    wall_color = [0,0,0] # Black
    left_color = [0,255,0] # Green
    right_color = [0,0,255] # Blue

    # Add shapes to the domain to create the corridor
    # Add the top horizontal wall to the domain
    line = Line2D( [0.0,width],[0.2*height,0.2*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom horizontal wall to the domain
    line = Line2D( [0.0,width],[0.8*height,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add a vertical exit on the left of the domain
    line = Line2D( [0.4*width,0.4*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=left_color, fill_color=left_color)
    # Add a vertical exit on the right of the domain
    line = Line2D( [0.6*width,0.6*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=right_color, fill_color=right_color)

    # Initalise the domain
    dom.build_domain()

    # Create destination objects towards each exits
    left = Destination(name='left', colors=[left_color],
                        excluded_colors=[wall_color])
    dom.add_destination(left)
    right = Destination(name='right', colors=[right_color],
                        excluded_colors=[wall_color])
    dom.add_destination(right)

    print("===> Domain: ",dom)

    # Define the spawning parameters for each group of people
    groups = [{"nb": 50, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [0.35,width/2-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "right"},

            {"nb": 50, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [width/2+0.35,width-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "left"}]

    # Plot the domain, distance to walls within the domain, and desired velocity
    # to each destination
    dom.plot(id=1, title="Domain", savefig=True,
            filename="domain.png")
    dom.plot_wall_dist(id=2, step=20,
        title="Distance to walls and its gradient",
        savefig=True, filename="room_wall_distance.png")
    dom.plot_desired_velocity('left',id=3, step=20,
        title="Distance to the left destination and desired velocity",
        savefig=True, filename="left_desired_velocity.png")
    dom.plot_desired_velocity('right',id=4, step=20,
        title="Distance to the right destination and desired velocity",
        savefig=True, filename="right_desired_velocity.png")

    return (dom, groups)

def domainInit2():
    width = 100.0
    height = 7.5

    # Create domain object
    dom = Domain(name = 'corridor',  pixel_size = 0.1, xmin=0, ymin=0,
                width=int(width*10), height=int(height*10))

    # Define wall and exit colours
    wall_color = [0,0,0] # Black
    left_color = [0,255,0] # Green
    right_color = [0,0,255] # Blue

    ## Add shapes to the domain to build a corridor
    # Add the top horizontal wall to the domain
    line = Line2D( [0.0,width],[0.2*height,0.2*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom horizontal wall to the domain
    line = Line2D( [0.0,width],[0.8*height,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add a vertical exit on the left of the domain
    line = Line2D( [0.4*width,0.4*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=left_color, fill_color=left_color)
    # Add a vertical exit on the right of the domain
    line = Line2D( [0.6*width,0.6*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=right_color, fill_color=right_color)

    # Initalise the domain
    dom.build_domain()

    # Create destination objects towards each exits
    left = Destination(name='left', colors=[left_color],
                        excluded_colors=[wall_color])
    dom.add_destination(left)
    right = Destination(name='right', colors=[right_color],
                        excluded_colors=[wall_color])
    dom.add_destination(right)

    print("===> Domain: ",dom)

    # Define the spawning parameters for each group of people
    groups = [{"nb": 100, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [0.35,width/2-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "right"},

            {"nb": 100, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [width/2+0.35,width-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "left"}]

    # Plot the domain, distance to walls within the domain, and desired velocity
    # to each destination
    dom.plot(id=1, title="Domain", savefig=True,
            filename="domain.png")
    dom.plot_wall_dist(id=2, step=20,
        title="Distance to walls and its gradient",
        savefig=True, filename="room_wall_distance.png")
    dom.plot_desired_velocity('left',id=3, step=20,
        title="Distance to the left destination and desired velocity",
        savefig=True, filename="left_desired_velocity.png")
    dom.plot_desired_velocity('right',id=4, step=20,
        title="Distance to the right destination and desired velocity",
        savefig=True, filename="right_desired_velocity.png")

    return (dom, groups)

def domainInit3():
    width = 100.0
    height = 7.5

    # Create domain object
    dom = Domain(name = 'corridor',  pixel_size = 0.1, xmin=0, ymin=0,
                width=int(width*10), height=int(height*10))

    # Define wall and exit colours
    wall_color = [0,0,0] # Black
    left_color = [0,255,0] # Green
    right_color = [0,0,255] # Blue

    ## Add shapes to the domain to build a corridor
    # Add the top horizontal wall to the domain
    line = Line2D( [0.0,width],[0.2*height,0.2*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom horizontal wall to the domain
    line = Line2D( [0.0,width],[0.8*height,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add a vertical exit on the left of the domain
    line = Line2D( [0.4*width,0.4*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=left_color, fill_color=left_color)
    # Add a vertical exit on the right of the domain
    line = Line2D( [0.6*width,0.6*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=right_color, fill_color=right_color)
    # Add the top vertical wall of the door to the domain
    line = Line2D( [0.5*width,0.5*width],[0.2*height,0.5*height-1], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom vertical wall of the door to the domain
    line = Line2D( [0.5*width,0.5*width],[0.5*height+1,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)

    # Initalise the domain
    dom.build_domain()

    # Create destination objects towards both exits
    left = Destination(name='left', colors=[left_color],
                        excluded_colors=[wall_color])
    dom.add_destination(left)
    right = Destination(name='right', colors=[right_color],
                        excluded_colors=[wall_color])
    dom.add_destination(right)

    print("===> Domain: ",dom)

    # Define the spawning parameters for each group of people
    groups = [{"nb": 25, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [0.35,width/2-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "right"},

            {"nb": 25, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [width/2+0.35,width-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "left"}]

    # Plot the domain, distance to walls within the domain, and desired velocity
    # to each destination
    dom.plot(id=1, title="Domain", savefig=True,
            filename="domain.png")
    dom.plot_wall_dist(id=2, step=20,
        title="Distance to walls and its gradient",
        savefig=True, filename="room_wall_distance.png")
    dom.plot_desired_velocity('left',id=3, step=20,
        title="Distance to the left destination and desired velocity",
        savefig=True, filename="left_desired_velocity.png")
    dom.plot_desired_velocity('right',id=4, step=20,
        title="Distance to the right destination and desired velocity",
        savefig=True, filename="right_desired_velocity.png")

    return (dom, groups)

def domainInit4():
    width = 100.0
    height = 7.5

    # Create domain object
    dom = Domain(name = 'corridor',  pixel_size = 0.1, xmin=0, ymin=0,
                width=int(width*10), height=int(height*10))

    # Define wall and exit colours
    wall_color = [0,0,0] # Black
    left_color = [0,255,0] # Green
    right_color = [0,0,255] # Blue

    ## Add shapes to the domain to build a corridor
    # Add the top horizontal wall to the domain
    line = Line2D( [0.0,width],[0.2*height,0.2*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom horizontal wall to the domain
    line = Line2D( [0.0,width],[0.8*height,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add a vertical exit on the left of the domain
    line = Line2D( [0.4*width,0.4*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=left_color, fill_color=left_color)
    # Add a vertical exit on the right of the domain
    line = Line2D( [0.6*width,0.6*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=right_color, fill_color=right_color)
    # Add the top vertical wall of the door to the domain
    line = Line2D( [0.5*width,0.5*width],[0.2*height,0.5*height-1], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom vertical wall of the door to the domain
    line = Line2D( [0.5*width,0.5*width],[0.5*height+1,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)

    # Initalise the domain
    dom.build_domain()

    # Create destination objects towards both exits
    left = Destination(name='left', colors=[left_color],
                        excluded_colors=[wall_color])
    dom.add_destination(left)
    right = Destination(name='right', colors=[right_color],
                        excluded_colors=[wall_color])
    dom.add_destination(right)

    print("===> Domain: ",dom)

    # Define the spawning parameters for each group of people
    groups = [{"nb": 50, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [0.35,width/2-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "right"},

            {"nb": 50, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [width/2+0.35,width-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "left"}]

    # Plot the domain, distance to walls within the domain, and desired velocity
    # to each destination
    dom.plot(id=1, title="Domain", savefig=True,
            filename="domain.png")
    dom.plot_wall_dist(id=2, step=20,
        title="Distance to walls and its gradient",
        savefig=True, filename="room_wall_distance.png")
    dom.plot_desired_velocity('left',id=3, step=20,
        title="Distance to the left destination and desired velocity",
        savefig=True, filename="left_desired_velocity.png")
    dom.plot_desired_velocity('right',id=4, step=20,
        title="Distance to the right destination and desired velocity",
        savefig=True, filename="right_desired_velocity.png")

    return (dom, groups)

def domainInit5():
    width = 100.0
    height = 7.5

    # Create domain object
    dom = Domain(name = 'corridor',  pixel_size = 0.1, xmin=0, ymin=0,
                width=int(width*10), height=int(height*10))

    # Define wall and exit colours
    wall_color = [0,0,0] # Black
    exit_color = [0,0,255] # Blue

    ## Add shapes to the domain to build a corridor
    # Add the top horizontal wall to the domain
    line = Line2D( [0.0,width],[0.2*height,0.2*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom horizontal wall to the domain
    line = Line2D( [0.0,width],[0.8*height,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add a vertical exit on the right of the domain
    line = Line2D( [0.6*width,0.6*width],[0.2*height,0.8*height], 2)
    dom.add_shape(line, outline_color=exit_color, fill_color=exit_color)
    # Add the top vertical wall of the door to the domain
    line = Line2D( [0.5*width,0.5*width],[0.2*height,0.5*height-0.6], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)
    # Add the bottom vertical wall of the door to the domain
    line = Line2D( [0.5*width,0.5*width],[0.5*height+0.6,0.8*height], 2)
    dom.add_shape(line, outline_color=wall_color, fill_color=wall_color)

    # Initalise the domain
    dom.build_domain()

    # Create destination objects towards the exit
    right = Destination(name='right', colors=[exit_color],
                        excluded_colors=[wall_color])
    dom.add_destination(right)

    print("===> Domain: ",dom)

    # Define the spawning parameters for the group of people
    groups = [{"nb": 100, "radius_distribution": ["uniform",0.25,0.35],
            "velocity_distribution":["normal",1.34,0.26],
            "box": [0.35,width/2-0.35,0.2*height+0.35,0.8*height-0.35],
            "destination": "right"}]

    # Plot the domain, distance to walls within the domain, and desired velocity
    # to the destination
    dom.plot(id=1, title="Domain", savefig=True,
            filename="domain.png")
    dom.plot_wall_dist(id=2, step=20,
        title="Distance to walls and its gradient",
        savefig=True, filename="room_wall_distance.png")
    dom.plot_desired_velocity('right',id=4, step=20,
        title="Distance to the right destination and desired velocity",
        savefig=True, filename="right_desired_velocity.png")

    return (dom, groups)

def simulation(domains, groups_set, seed=1):
    # Define parameters
    dmax = 0.1
    F = 2000.0
    Fwall = 2000.0
    lambda_ = 0.5
    delta = 0.08
    kappa = 120000.0
    eta = 240000.0
    tau = 0.5
    mass = 80.0
    Tf = 90.0 # Final time
    t = 0.0 # Current time
    dt = 0.005 # Timestep

    counter = 0 # Plot counter
    cc = 0 # Simulation loops since plot counter
    drawper = 200.0 # Simulation loops required for plot
    draw = True

    print("===> Seed = "+str(seed))

    all_people = {}
    for name, groups in groups_set.items():
        dom = domains[name]

        # Initialise all people within the domain
        people = people_initialization(dom, groups, dt, dmin_people=0.0,
                                        dmin_walls=0.0, seed=seed,
                                        itermax=10, verbose=True,
                                        projection_method="cvxopt")
        I, J, Vd = dom.people_desired_velocity(people["xyrv"],
                                                people["destinations"])
        people["Vd"] = Vd
        for ip,pid in enumerate(people["id"]):
            people["paths"][pid] = people["xyrv"][ip,:2]
        contacts = None
        all_people[name] = people

    print("===> All people = ",all_people)

    # Create the figure so that the text size is outputted correctly for
    # future plots
    plot_people(5, dom, [], [],[],time=0.0)
    plt.pause(0.01)

    while (t<Tf):
        print("\n===> Time = "+str(t))

        for name, dom in domains.items():
            people = all_people[name]
            print("===> Compute desired velocity for the domain ", name)
            I, J, Vd = dom.people_desired_velocity(people["xyrv"],
                                                    people["destinations"])
            people["Vd"] = Vd
            people["I"] = I
            people["J"] = J

        # Duplicate people who exist in transit boxes to their destination
        # domain (not used with single domains)
        print("===> Find people who have to be duplicated")
        virtual_people = find_duplicate_people(all_people, domains)
        print("     virtual_people : ",virtual_people)

        for name, dom in domains.items():
            people = all_people[name]
            try:
                xyrv = np.concatenate((people["xyrv"],
                    virtual_people[name]["xyrv"]))
                Vd = np.concatenate((people["Vd"],
                    virtual_people[name]["Vd"]))
                Uold = np.concatenate((people["Uold"],
                    virtual_people[name]["Uold"]))
            except:
                xyrv = people["xyrv"]
                Vd = people["Vd"]
                Uold = people["Uold"]

            print("===> Compute social forces for domain ", name)

            # If there are people in the domain
            if (xyrv.shape[0]>0):

                # If there are multiple people occupying the same position
                if (np.unique(xyrv, axis=0).shape[0] != xyrv.shape[0]):
                    print("===> ERROR : There are two identical lines in the")
                    print("             array xyrv used to determine the \
                        contacts between")
                    print("             individuals and this is not normal.")
                    sys.exit()

                contacts = compute_contacts(dom, xyrv, dmax)
                print("     Number of contacts: ",contacts.shape[0])
                Forces = compute_forces( F, Fwall, xyrv, contacts, Uold, Vd,
                                         lambda_, delta, kappa, eta)
                nn = people["xyrv"].shape[0]
                all_people[name]["U"] = dt*(Vd[:nn,:]-Uold[:nn,:])/tau + \
                                              Uold[:nn,:] + \
                                              dt*Forces[:nn,:]/mass
                # Calculate the velocity of virtual people for plotting purposes
                # (this velocity is not used for any future computation)
                virtual_people[name]["U"] = dt*(Vd[nn:,:]-Uold[nn:,:])/tau + \
                                                  Uold[nn:,:] + \
                                                  dt*Forces[nn:,:]/mass

                all_people[name], _ = move_people(t, dt, all_people[name], {})
            else:
                nn = 0

            if(draw):
                counter += 1
                # Colour people going towards the left exit in a different
                # colour
                colors =  np.zeros(nn)
                for i in range(nn):
                    if people["destinations"][i] == "left":
                        colors[i] = 1


                plot_people(5, dom, all_people[name], contacts,
                            colors, virtual_people=virtual_people[name], time=t,
                            plot_people=True, plot_contacts=True,
                            plot_paths=False, plot_velocities=False,
                            plot_desired_velocities=False, plot_sensors=False,
                            sensors={}, savefig=True,
                            filename=name+"_fig_"+str(counter)+'.png')
                plt.pause(0.01)

                if(t==0):
                    plot_people(6, dom, all_people[name], contacts,
                                colors, virtual_people=virtual_people[name],
                                time=t, plot_people=True, plot_contacts=True,
                                plot_paths=False, plot_velocities=False,
                                plot_desired_velocities=True,
                                plot_sensors=False, sensors={}, savefig=True,
                                filename='initial_velocities.png')
                    plt.pause(0.01)


            # Update people destinations
            all_people = people_update_destination(all_people,domains,
                                                    dom.pixel_size)
            # Update previous velocities
            all_people[name]["Uold"] = all_people[name]["U"]

            # Print the number of persons for each domain
            print("===> Domain " + name + " nb of persons = ",
                    all_people[name]["xyrv"].shape[0])

            t += dt
            cc += 1
            if (cc>=drawper):
                draw = True
                cc = 0
            else:
                draw = False

    plt.ioff()
    plt.show()
    sys.exit()

def main(i):
    if(i==1):
        dom, groups = domainInit()
        groups_set = {}
        groups_set["corridor"] = groups
        simulation({"corridor":dom},groups_set,2)
    if(i==2):
        dom, groups = domainInit2()
        groups_set = {}
        groups_set["corridor"] = groups
        simulation({"corridor":dom},groups_set)
    if(i==3):
        dom, groups = domainInit3()
        groups_set = {}
        groups_set["corridor"] = groups
        simulation({"corridor":dom},groups_set)
    if(i==4):
        dom, groups = domainInit4()
        groups_set = {}
        groups_set["corridor"] = groups
        simulation({"corridor":dom},groups_set)
    if(i==5):
        dom, groups = domainInit5()
        groups_set = {}
        groups_set["corridor"] = groups
        simulation({"corridor":dom},groups_set)

main(5)
