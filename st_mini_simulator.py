import streamlit as st
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import io 
import matplotlib
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import csv

st.title("Mini Simulator for Reaction Kinetics, Heat Transfer and Greenhouse Gas Impact ")
type_of_cal = st.selectbox("Heat Transfer, Reaction Kinetics or Greenhouse Gas Impact?:", ["---Please select---", "Heat Transfer", "Reaction Kinetics", "Greenhouse Gas Impact"])

def type_of_calculation():
    if type_of_cal == "---Please select---":
        st.subheader("Hello! This is a lightweight, interactive tool to simulate selected introductory concepts in Chemical Engineering.")
        st.markdown("#### Here, you can explore: \n ##### 1. Reaction kinetics (supports graphing for Zero, First and Second-order reactions; custom molar ratios.) \n ##### 2. Heat Transfer Calculations (supports graphing to compare materials in phase changes, conduction and radiation) \n ##### 3. Greenhouse Gas Impact (calculates CO₂ and methane emissions, supports graphing to compare properties of fuels)")
        st.markdown("#### Enjoy!")

    if type_of_cal == "Reaction Kinetics":
        st.subheader("Calculating the Rate Constant")
        def get_k():
            calculate_k = st.selectbox("How would you like to calculate k? (Rate and concentration/Arrhenius equation)?", ["---Please select---", "Rate and concentration", "Arrhenius equation"])
            if calculate_k == "Arrhenius equation":
                A = st.number_input("A (frequency factor):", min_value=0.0)
                Ea = st.number_input("Ea (activation energy in J/mol):", min_value=0.0)
                T = st.number_input("Temperature (K):", min_value = 0.0)
                R = 8.31
                k = None
                if st.button("Calculate rate constant"):
                    if T == 0:
                        st.warning("Temperature cannot be zero!")
                    else:
                        k = A * math.exp(-Ea / (R * T))
                        st.session_state.k = k
                        st.success(f"Rate constant = {st.session_state.k:.6f}")
                if "k" in st.session_state:
                    st.info(f"Rate constant = {st.session_state.k:.6f}")
                return k 
            if calculate_k == "Rate and concentration":
                k = None
                rate_str = st.text_input("Rate in molL⁻¹s⁻¹: ")
                try:
                    rate = float(rate_str)
                except ValueError:
                    st.warning("Please enter a valid number! ")
                    return
                concentration = st.number_input("Concentration: ")
                if rate <= 0 or concentration <=0:
                    st.warning("Rate and concentration cannot be less than 0!")
                    return
                order_selected = st.selectbox("Order with respect to the substance? Please select from the following:", ["---Please select---", 0.0, 1.0, 2.0])
                if order_selected == "---Please select---":
                    return
                k = rate/(concentration ** order_selected)
                st.session_state.k = k

                if "k" in st.session_state:
                    st.info(f"Rate constant = {st.session_state.k:.6f}")
                return k

        def get_conc():
            st.subheader("Information for Plotting Graph")
            range_or_data = st.selectbox("Plot graph using manually inputted data or auto-generated time values based on the range given?", ["---Please select---", "Manually inputted data", "Auto-generate time values based on range given"])
            k = st.session_state.get("k", None)
            A0 = st.number_input("Initial concentration (mol/L):", min_value=0.0)
            order = st.selectbox("Order of Reaction:", ["---Please select---", "Zero", "First", "Second"])
            if order == "---Please select---":
                st.warning("Please select an order!")
                return
            else:
                st.info(f"Order selected: {order} Order")
            st.markdown("##### Product's Information")
            ratio = st.number_input("What is the ratio of the product to the reactant?", min_value = 0.0)

            if range_or_data == "Manually inputted data":
                st.markdown("##### Fill in Time Values to Calculate Concentrations")
                n = st.number_input("Number of times to calculate concentration:", min_value=1.0)
                time  = []
                R_conc = []
                P_conc = []
                if k is not None:
                    for i in range(int(n)):
                        t = st.number_input(f"Time: {i+1}", min_value=0.0, key=f"time_{i}")
                        if order == "Second":
                            if A0 == 0:
                                st.warning("Initial Concentration cannot be zero! ")
                                break
                            else:
                                inverse_R = 1/A0 + k*t
                                R_concentration = 1/inverse_R

                        elif order == "First":
                            R_concentration = A0*math.exp(-k*t)

                        if order == "Zero":
                            R_concentration = A0 - k*t
                            if R_concentration < 0:
                                st.warning("Reactant's concentration cannot be less than 0!")
                                continue
                        time.append(round(t, 3))

                        P_concentration = (A0 - R_concentration)*ratio
                        st.success(f"[Reactant] = {R_concentration:.3f}mol/L")
                        st.success(f"[Product] = {P_concentration:.3f}mol/L")
                        R_conc.append(round(R_concentration,3))
                        P_conc.append(round(P_concentration, 3))
                else:
                    st.warning("Please calculate rate constant first!")
                if order == "Second":
                    st.write("The half-life of second order reactions increases as the reactant's concentration decreases.")
                elif order == "First":
                    if k is not None:
                        half_life = math.log(2)/k
                        st.write(f"Since the half-life of first-order reactions is constant, half-life = {half_life}s")
                    else:
                        return
                elif order == "Zero":
                    st.write("The half-life of zero order reactions decreases as the reactant's concentration decreases.")

                st.subheader("Time-Concentration Data (Manual Input)")
                st.write(
                pd.DataFrame(
                    {   
                        "time (s)": (time),
                        "Concentration of Reactant (mol/L)": (R_conc),
                        "Concentration of Product (mol/L)": (P_conc),
                    }
                )
            )
            st.subheader("Concentration Vs. Time Graph")
            #conc and time cannot be negative!
            if range_or_data == "Auto-generate time values based on range given":
                if k is None:
                    st.warning("Please calculate rate constant first!")
                    return
                fig,ax = plt.subplots(figsize = (12,6))
                time_min = st.number_input("Minimum time:", min_value = 0.0)
                time_max = st.number_input("Maximum time", min_value = 0.0)
                num_points = 100
                times = np.linspace(time_min, time_max, num_points)
                plt.figure(figsize=(12,6))
                if order == "Second":
                    inverse_R = 1/A0 + k*times
                    R_concentration = 1/inverse_R
                elif order == "First":
                    R_concentration = A0*np.exp(-k*times)
                elif order == "Zero":
                    R_concentration = A0 - k*times
                P_concentration = (A0 - R_concentration)*ratio 

                ax.plot(times, R_concentration, label ="Reactant")
                ax.plot(times, P_concentration, label = "Product")
                ax.set_xlabel("time(s)")
                ax.set_ylabel("Concentration(mol/L)")
                plot_title = st.text_input("Graph's title: ")
                ax.set_title(f"{plot_title}")
                ax.grid(True)
                ax.legend()
                st.pyplot(fig)

            if range_or_data =="Manually inputted data":
                style = st.selectbox("Animated or static?", ["---Please select---", "Static","Animated"])
                if style == "Static":
                    fig, ax = plt.subplots()
                    y = np.array(R_conc)

                    ax.plot(time, y, marker='o', linestyle='-', color = "r", label='reactant')
                    ax.plot(time, P_conc, marker='*', linestyle="-", color = "b", label='product')
                    ax.set_xlabel("time(s)")
                    ax.set_ylabel("Concentration(mol/L)")
                    plot_title = st.text_input("Graph's title: ")
                    ax.set_title(f"{plot_title}")
                    ax.grid(True)
                    ax.legend()
                    st.pyplot(fig)

                    st.subheader("Download File")
                    format = st.selectbox("File Format:", ["---Please select---", "png", "jpeg", "pdf"])
                    if format != "---Please select---":
                        buf = io.BytesIO()
                        fig.savefig(buf, format = format)
                        buf.seek(0)
                    else:
                        st.warning("Please select a file format.")

                    file_name = st.text_input("File name: ", value = f"kinetics_plot.{format}")
                    if format == "png" or format == "jpeg":
                        st.download_button(label = "Download image", data = buf, file_name = f"{file_name}", mime=f"image/{format}")
                    elif format == "pdf":
                        st.download_button(label = "Download image", data = buf, file_name = f"{file_name}", mime="application/pdf")
        
                elif style == "Animated":
                    if not time:
                        return
                    else:
                        fig, ax = plt.subplots()

                        x_min = min(time)
                        x_max = max(time)
                        y_min = min(min(R_conc), min(P_conc))
                        y_max = max(max(R_conc), max(P_conc))
                        ax.set_xlim(x_min, x_max)
                        ax.set_ylim(y_min, y_max)
                        line1, = ax.plot([], [], lw = 2, marker='o', linestyle='-', color = "r", label='reactant')
                        line2, = ax.plot([], [], lw = 2, marker='*', linestyle="-", color = "b", label='product')
                        ax.set_xlabel("time(s)")
                        ax.set_ylabel("Concentration(mol/L)")
                        plot_title = st.text_input("Graph's title: ")
                        ax.set_title(f"{plot_title}")
                        ax.grid(True)
                        ax.legend()

                        def init(): 
                            line1.set_data([], []) 
                            line2.set_data([], [])
                            return line1, line2

                        def animate(i): 
                            x = time[:i+1] 
                            y1 = P_conc[:i+1]
                            y2 = R_conc[:i+1]
                            line1.set_data(x, y1) 
                            line2.set_data(x, y2)
                            return line1, line2

                        anim = FuncAnimation(fig, animate, 
                                    init_func = init, 
                                    frames = len(time), 
                                    interval = 1000, 
                                    blit = True) 

                        anim.save("Kinetics_animation.gif", 
                        writer = "pillow", fps = 2)
                        st.image("Kinetics_animation.gif", caption="Kinetics Animation!")
                else:
                    return
        def helper():
            get_k()
            get_conc()
        helper()
    if type_of_cal == "Heat Transfer":
        st.markdown("##### Types of calculations supported: ")
        st.info("Note: The item in brackets shows the parameter this option will calculate.")
        st.info("Note: 📈 means graphing is available for this option.")
        table_type = ["1. Sensible Heating or Cooling (heat removed/added or final temperature)",
             "2. Phase Changes with known heat of vaporisation (Rate of heat transfer)",
             "3. Phase Changes with no known heat of vaporisation 📈 (Heat of vaporisation + Rate of heat transfer)",
             "4. Chemical Reactions (heat added/removed)",
             "5. Heat Transfer Area of Exchanger (area of exchanger)",
             "6. Conduction - Fourier’s law 📈 (rate of heat transfer)",
             "7. Radiation - Stefan-Boltzmann Law 📈 (rate of heat transfer)",
             "8. Convective Heat Transfer (rate of heat transfer)"]
        table_ = []
        for eqn in table_type:
            table_.append({"Type of Calculation":eqn})
        df = pd.DataFrame(table_) 
        st.dataframe(df)
        
        calculation = st.selectbox("Choose the type of calculation from the following:", ["---Please select---", "1. Sensible Heating or Cooling",
             "2. Phase Changes with known heat of vaporisation",
             "3. Phase Changes with no known heat of vaporisation",
             "4. Chemical Reactions",
             "5. Heat Transfer Area of Exchanger",
             "6. Conduction - Fourier’s law",
             "7. Radiation - Stefan-Boltzmann Law",
             "8. Convective Heat Transfer"])

        if calculation == "1. Sensible Heating or Cooling":
            n = st.number_input("Number of feed streams: ", min_value = 1.0) #no unit
            hc_product = st.number_input("Specific Heat capacity of product stream in kJ/(kg K) : ", min_value = 1e-6)
            if hc_product == 1e-6:
                st.warning("Specific heat capacity of product stream cannot be empty!") 
                return
            T_ref = 25
            ṁ_list = []
            feed_stream_list = []
            st.subheader("Calculate Heat Flow or Product Temperature")
            type_of_cal_SH = st.selectbox("(1) Calculate amount of heat removed/added or (2) temperature of product?", ["---Please select---", "Amount of heat removed/added", "Temperature of product"])
            for i in range(int(n)): #asking for information on the feed stream
                st.markdown(f"##### Information for feed stream {i + 1}")
                ṁ = st.number_input(f"Mass flow rate of feed stream {i + 1} in kg/unit time: ") #kg/min or kg/unit time #must not be in moles!
                ṁ_list.append(ṁ)
                T = st.number_input(f"Temperature of feed stream {i + 1} in celsius: ") #unit: celsius
                hc = st.number_input(f"Specific heat capacity of feed stream {i + 1} in kJ/(kg K): ") #unit: kJ/(kg K)
                sum_feed_stream = hc * ṁ * (T - T_ref)
                feed_stream_list.append(sum_feed_stream)
            sum_ṁ = sum(ṁ_list) #Mass flow rate of product
            st.info(f"Steady-state Total Mass Balance (Mass flow rate of product): {sum_ṁ} kg/unit time")
            if ṁ == 0:
                st.warning("Mass flow rate of feed stream cannot be empty!")
                return
            if T == 0:
                st.warning("Temperature of feed stream cannot be empty!")
                return
            if hc == 0:
                st.warning("Specific heat capacity of feed stream cannot be empty!")
                return
            if type_of_cal_SH == "Amount of heat removed/added": #Calculate amount of heat removed/added
                T_product = st.number_input("Temperature of product stream in celsius: ") #unit: celsius
                if T_product == 1e-6:
                    st.warning("Temperature of product stream cannot be empty ")
                    return
                Q̇ = sum_ṁ * hc_product * (T_product - T_ref) - sum(feed_stream_list)
                if "-" in str(Q̇):
                    st.info(f"Heat removed: {round(Q̇, 2)} kJ/unit time") #unit: depends on the unit of time used in mass flow rate
                    return
                if Q̇ == 0:
                    st.info("It is an adiabatic process as Q̇ = 0!")
                    return
                else:
                    st.info(f"Heat added: {round(Q̇, 2)} kJ/unit time")
            elif type_of_cal_SH == "Temperature of product": #temperature of product
                if sum_ṁ == 0:
                    st.warning("Please calculate the Mass flow rate of product first!")
                    return
                if hc_product == 0:
                    st.warning("The specific heat capacity of the product stream cannot be zero!")
                    return
                T_product = sum(feed_stream_list)/(sum_ṁ * hc_product) + T_ref
                st.info(f"Temperature of product: {round(T_product, 2)} celsius")
                #cannot have curly braces here around T_product! It becomes a set

        if calculation == "2. Phase Changes with known heat of vaporisation":
            m = st.number_input("Mass flow rate (in kg/min): ", min_value = 1e-6)
            h_vap = st.number_input("Heat of vaporisation (in J/kg): ") #unit: in J/kg
            result = m * h_vap
            st.info(f"Q̇ (Rate of transfer of energy) = {result} J/min")

        if calculation == "3. Phase Changes with no known heat of vaporisation":
            ask_user = st.selectbox("(1) Temperature against Pressure graph using the Clausius-Clapeyron Equation or (2) Calculate the heat of vaporisation and rate of heat transfer?", 
                                    ["---Please select---", "Temperature against Pressure graph using the Clausius-Clapeyron Equation", "Calculate the heat of vaporisation and rate of heat transfer"])
            if ask_user == "Temperature against Pressure graph using the Clausius-Clapeyron Equation":
                graph_type = st.selectbox("(1) Explore common substances or (2) plot your own graph?", ["---Please select---","Explore common substances", "Plot your own graph"], key = "select_type")
                R = 8.31 #Molar gas constant #in J/(K mol)

                if graph_type == "Explore common substances":
                    table = [] #preparing table for appending
                    with open("list_of_boiling_points_and_heat_of_vap.csv", "r", encoding = "utf-8") as file: #open the csv file
                        reader = csv.DictReader(file)
                        for row in reader: #arranging the components inside the csv into a table
                            table.append({"Compound": row["Compound"],
                                            "Boiling Point at Normal Pressure(K)":row["Boiling Point at Normal Pressure(K)"],
                                            "Heat of Vaporization(J/mol)": row["Heat of Vaporization(J/mol)"]})

                    st.subheader("Boiling Point and Heat of Vaporisation Table (Standard Pressure)")
                    df = pd.DataFrame(table) 
                    st.dataframe(df)

                    table_substance = {}#preparing a dict for putting what the user wants into a dict
                    compound_names = [name["Compound"] for name in table] #if separate for loop and use indentation instead - over-writing compound names for each iteration
                    common_substance = st.multiselect("Choose a material:", compound_names) #returns a list
                    st.subheader("Graphs")
                    if not common_substance:
                        st.warning("To display the graph, please select at least one item from the table!")
                        return
                    for substance in common_substance:
                        for item in table:
                            if substance == item["Compound"]:
                                #must add float as substance["..."] is a string!
                                T1 = float(item["Boiling Point at Normal Pressure(K)"]) #in Kelvin
                                L = float(item["Heat of Vaporization(J/mol)"]) #in J/mol
                                P1 = 101325 #in pa
                                #cannot use "append" bcz it is a dict. Use: dict["key"] = value
                                table_substance[item["Compound"]] = {"P1":P1,
                                                                            "T1":T1,
                                                                            "L": L} #all substances are now in the dict!
                    fig1, ax1 = plt.subplots(figsize = (12,6))
                    #for the curve
                    T_min = 250 #in K
                    T_max = 400 #in K
                    num_points = 100
                    #linspace() function in NumPy returns an array of evenly spaced numbers over a specified range
                    temperatures = np.linspace(T_min, T_max, num_points)
                    T2 = temperatures
                    #to access the items in the dictionary(table_substance)
                    for substance_s, params in table_substance.items():
                        P1 = params["P1"]
                        T1 = params["T1"]
                        L = params["L"]
                        #for the first graph (curve)
                        P2 = P1 * np.exp(-L/R * (1/T2 - 1/T1))
                        ax1.plot(temperatures, P2, label=f"{substance_s}")
                    #curve
                    ax1.set_xlabel("Temperature(K)")
                    ax1.set_ylabel("Pressure(Pa)")
                    ax1.set_title("Pressure(Pa) VS Temperature(K)")
                    ax1.grid(True)
                    ax1.legend()
                    st.pyplot(fig1)

                    fig2, ax2 = plt.subplots(figsize = (12,8))
                    T_min_linear = 200 #in K
                    T_max_linear = 300 #in K
                    num_points_linear = 200
                    temperatures_linear = np.linspace(T_min_linear, T_max_linear, num_points_linear)
                    T2_linear = temperatures_linear

                    for substance_s, params in table_substance.items():
                        P1 = params["P1"]
                        T1 = params["T1"]
                        L = params["L"]
                        #for the first graph (curve)
                        P2_linear = P1 * np.exp(-L/R * (1/T2_linear - 1/T1))
                        ax2.plot(1/T2_linear, np.log(P2_linear), label=f"{substance_s}")
                        #plot the graph

                    ax2.set_xlabel("1/Temperature(K)")
                    ax2.set_ylabel("ln(Pressure(Pa))")
                    ax2.set_title("ln(Pressure(Pa)) VS 1/Temperature(K)")
                    ax2.grid(True)
                    ax2.legend()
                    st.pyplot(fig2)
                    #input file name

                #use real experimental data
                if graph_type == "Plot your own graph":
                    temp_graph_type_2 = []
                    vp_graph_type_2 = [] #appending user's data into lists
                    st.subheader("Data Collection for Graphing")
                    number = st.number_input("How many pieces of data? ", min_value=1.0, key = "s")
                    for i in range(int(number)):
                        temp_2 = st.number_input("Temperature in K: ", min_value=1e-6, key = f"t{i}")
                        temp_graph_type_2.append(1/temp_2)
                        vapor_pressure = st.number_input("Vapour pressure: ", min_value=1e-6, key=f"g{i}")
                        vp_graph_type_2.append(math.log(vapor_pressure))
                    #plots graph with line of best fit
                    fig, ax = plt.subplots()
                    y = np.array(vp_graph_type_2)
                    x = np.array(temp_graph_type_2)
                    slope, intercept = np.polyfit(temp_graph_type_2, y, 1)
                    #np.array() is for number-crunching, allows proper multiplication (or maths in general)
                    #normal lists [1, 2, 3] * 2 = [1, 2, 3, 1, 2, 3]
                    line = slope * x + intercept
                    ax.scatter(x, y, color="blue", marker = "*")
                    st.markdown("##### Graph's Information")
                    label_name = st.text_input("Enter a label for the line of best fit:")
                    ax.plot(x, line, color="red", label=f"Line of Best Fit, {label_name}")
                    #graph's info
                    ax.set_xlabel("1/Temperature(K)")

                    unit = st.text_input("Unit of Pressure: ")
                    ax.set_ylabel(f"ln(Vapour Pressure({unit}))")
                    plot_title = st.text_input("Graph's title: ")
                    st.subheader("Graph")
                    ax.set_title(f"{plot_title}")
                    ax.grid(True)
                    #no need plt.legend() if there is no label
                    ax.legend()
                    st.pyplot(fig)

            if ask_user == "Calculate the heat of vaporisation and rate of heat transfer":
                st.markdown("##### Temperature")
                T1 = st.number_input("Temperature (T1) in K: ", min_value = 0) #in Kelvin
                T2 = st.number_input("Temperature (T2) in K: ", min_value = 0) #in Kelvin
                if T1 == 0 or T2 == 0:
                    st.warning("Temperature cannot be less than or equal to 0!")
                    return
                st.markdown("##### Vapour Pressure")
                P1 = st.number_input(f"Vapour pressure at {T1}K: ", min_value = 0) #any unit
                P2 = st.number_input(f"Vapour pressure at {T2}K: ", min_value = 0) #any unit
                if P1 == 0 or P2 == 0:
                    st.warning("Vapour pressure cannot be less than or equal to zero!")
                    return
                st.markdown("##### Mass Flow Rate")
                m = st.number_input("Mass flow rate in mol/unit time: ", min_value = 0) # unit: mol/unit time
                R = 8.31 #in J/(K mol)
                L = (R*math.log(P1/P2))/(1/T2 - 1/T1) #calculate Heat of vaporisation
                st.info(f"Heat of vaporisation: {L} J/mol")
                Q̇ = m * L #calculate Rate of transfer of energy
                st.info(f"Q̇ (Rate of transfer of energy): {round(Q̇/1000,2)} kJ/unit time")

        if calculation == "4. Chemical Reactions":
            percentage = st.number_input("Percentage of reactant converted to product in %: ", min_value=1e-6)
            r = st.number_input("Rate at which limiting reactant is consumed (in mol/unit time): ", min_value=1e-6) #unit: in mol/unit time
            enthalpy_change = st.number_input("ΔH (Heat of reaction) in J/mol: ")
            if enthalpy_change == 0:
                return
            if percentage == 1e-6 or r == 1e-6  or enthalpy_change == 1e-6:
                return
            if percentage > 100:
                st.warning("The percentage cannot be greater than 100!")
                return
            Q̇ = (percentage/100) *r* enthalpy_change
            if "-" in str(Q̇):
                st.info(f"Heat removed: {round(Q̇/1000,3)} kJ/unit time. Hence, since the reaction is exothermic, {round(Q̇/1000,3)} kJ/unit time needs to be removed to keep the temperature constant in the reactor.") #unit: depends on the unit of time used in mass flow rate
                return
            else:
                st.info(f"Heat added: {round(Q̇/1000,3)} kJ/unit time. Hence, since the reaction is endothermic, {round(Q̇/1000, 3)} kJ/unit time needs to be added to keep the temperature constant in the reactor.")

        if calculation == "5. Heat Transfer Area of Exchanger":
            Uo_dict = {"Saturated vapor: Boiling liquid": 250, #unit: Btu/hr ft² °F
               "Saturated vapor: Flowing liquid": 150,
               "Saturated vapor: Vapor": 20,
               "Liquid: Liquid": 50,
               "Liquid: Gas or Gas:Liquid": 20,
               "Gas: Gas": 10,
               "Vapor(partial condenser): Liquid":30}

            st.subheader("Table: Approximate Values of Uo")
            table_Uo = []
            for item in Uo_dict:
                table_Uo.append({"Hot Stream: Cold Stream": item, "Uo, Btu/hr ft² °F": Uo_dict[item]})

            df = pd.DataFrame(table_Uo)
            st.dataframe(df)
            Uo = st.number_input("Based on the table above or your own values, overall heat transfer coefficient (Uo, in Btu/hr ft² °F): ") #unit: Btu/hr ft² °F
            Q_duty = st.number_input("Rate of heat transfer in Btu/min: ") #Unit: Btu/min #can be -ve, when heat is transferred out of a system
            delta_T1 = st.number_input("ΔT1: ") #unit:°F
            delta_T2 = st.number_input("ΔT2: " ) #unit:°F
            if delta_T1 < 0 and delta_T2 > 0 or delta_T2 < 0 and delta_T1 > 0:
                st.warning("ΔT1 and ΔT2 must not have opposite signs! Please try again." )
                return
            if Uo == 0:
                st.warning("Overall heat transfer coefficient cannot be zero!")
                return
            if delta_T2 == 0:
                st.warning("ΔT2 canot be zero!")
                return
            T_lm = (delta_T1 - delta_T2)/math.log(delta_T1/delta_T2)
            print(f"ΔTlog mean = {T_lm} °F")
            A = (Q_duty/(Uo * T_lm))*60
            st.info(f"Area = {A} ft²") #Unit: ft²
            if A < 0:
                st.warning("The area cannot be less than zero!")

        if calculation == "6. Conduction - Fourier’s law":
            ask_user = st.selectbox("(1) Temperature against Distance graph of different materials or (2) calculate rate of energy transfer?", ["---Please select---", "Temperature against Distance graph of different materials", "Calculate rate of energy transfer"])
            if ask_user == "Calculate rate of energy transfer":
                st.subheader("Calculation using the Fourier's Law")
                T_high = st.number_input("Temperature of high-temperature region in K:", min_value = 1e-6) #unit: K
                T_low = st.number_input("Temperature of low-temperature region in K:", min_value=1e-6) #unit: K
                if T_high < 1e-6 or T_low < 1e-6:
                    st.warning("Temperatures cannot be less than 0 Kelvin!") 
                    return
                d = st.number_input("Thickness in m: ")#unit: m
                K_tc = st.number_input("Thermal conductivity in W/(m K): ") #unit: W/(m K)
                area_conduction = st.number_input("Cross-sectional area in m²: ") #unit: m²
                if d <= 0 or K_tc <= 0 or area_conduction <= 0:
                    st.warning("Thickness, Thermal conductivity and cross-sectional area cannot be <= 0!")
                    return
                rate_of_heat_transfer  = (K_tc * area_conduction * (T_high - T_low))/d
                st.info(f"Rate of energy transfer = {rate_of_heat_transfer} W")

            if ask_user == "Temperature against Distance graph of different materials":
                table = [] #preparing a list for appending
                with open("list_of_thermal_conductivities.csv","r", encoding = "utf-8") as file: #open csv
                    reader = csv.DictReader(file)
                    for row in reader:
                        table.append({"Material": row["Material"],
                                        "Thermal Conductivity (W·m⁻¹·K⁻¹) at atmospheric pressure and around 293 K (20 °C)":row["Thermal Conductivity (W·m⁻¹·K⁻¹) at atmospheric pressure and around 293 K (20 °C)"]})
                st.subheader("Table: Material and Thermal Conductivity")
                df = pd.DataFrame(table)
                st.dataframe(df)

                table_substance = {}
                material_names = (name["Material"] for name in table)
                common_substance = st.multiselect(f"Please select one or more items from the list below:", material_names) #returns list
                st.subheader("Graph")
                if not common_substance:
                    st.warning("To display the graph, please select at least one item from the table!")
                    return
                for substance in common_substance: #to ensure substances e.g. Aluminium and Aluminium nitride don't get mixed up
                    for item in table:
                        if substance == item["Material"]:
                            K = float(item["Thermal Conductivity (W·m⁻¹·K⁻¹) at atmospheric pressure and around 293 K (20 °C)"])
                            T_high = 75 #in celsius
                            area_conduction = 0.5 #in m^2
                            Q = 500 #rate #in W
                            table_substance[item["Material"]] = {"K":K}

                fig, ax = plt.subplots(figsize = (12,6))
                d_min = 0 #min thickness #in m
                d_max = 1.5 #max thickness in m
                num_points = 100
                distances = np.linspace(d_min, d_max, num_points) #returns an array of evenly spaced numbers over a specified range

                for substance_s, thermal_conductivity in table_substance.items():
                    K = thermal_conductivity["K"]
                    T_low = T_high - (Q / (K * area_conduction)) * distances # calculate T_low based on array of numbers
                    T_low = np.maximum(T_low, -273) #make sure lowest T doesn't go below -273 celsius (lowest possible temp)
                    ax.plot(distances, T_low, label=f"{substance_s}, K = {thermal_conductivity['K']} W·m⁻¹·K⁻¹") #plotting starts
                
                #graph's info
                ax.set_xlabel("Distance (m)")
                ax.set_ylabel("Temperature (°C)")
                plot_title = st.text_input("Graph's title: ")
                ax.set_title(f"{plot_title}")
                ax.grid(True)
                ax.legend()
                st.pyplot(fig)

        if calculation == "7. Radiation - Stefan-Boltzmann Law":
            ask_user = st.selectbox("(1) Power against Temperature graph of different materials or (2) calculate rate of energy transfer?", ["---Please select---", "Power against Temperature graph of different materials", "Calculate rate of energy transfer"])
            if ask_user == "Calculate rate of energy transfer":
                st.subheader("Calculation using the Stefan-Boltzmann Law")
                stefan_boltzmanns_constant = 5.67e-8 #unit: W/m²K⁴
                emissivity = st.number_input("Emissivity:") #no unit
                if emissivity > 1 or emissivity < 0:
                    st.warning("Emissivity can only be between 1 and 0!")
                    return
                surface_area = st.number_input("Surface Area in m²:", min_value = 1e-6) #unit: m²
                temperature = st.number_input("Temperature in K:", min_value=1e-6) #unit: K
                if surface_area <= 1e-6 or temperature <= 1e-6:
                    st.warning("Temperature and surface area cannot be <= 0!")
                    return
                heat_flow_rate = emissivity * stefan_boltzmanns_constant * surface_area * temperature ** 4
                st.info(f"Rate of energy transfer = {heat_flow_rate} W")

            if ask_user == "Power against Temperature graph of different materials":
                table = [] #preparing a list for appending
                with open("list_of_emissivities.csv","r") as file: #open csv
                    reader = csv.DictReader(file)
                    for row in reader:
                        table.append({"Material": row["Material"],
                                        "Emissivity at 300 K":row["Emissivity at 300 K"]})
                
                st.subheader("Table: Material and Emissivity")
                df = pd.DataFrame(table)
                st.dataframe(df)

                table_substance = {}
                material_names = (name["Material"] for name in table)
                common_substance = st.multiselect(f"Please select one or more items from the list below:", material_names) #returns list
                st.subheader("Graph")
                if not common_substance:
                    st.warning("To display the graph, please select at least one item from the table!")
                    return
                table_substance = {}
                for substance in common_substance:
                    for item in table:
                        if substance in item["Material"]:
                            emissivity = float(item["Emissivity at 300 K"])
                            stefan_boltzmanns_constant = 5.67e-8
                            surface_area = 1.0 #in m^2
                            table_substance[item["Material"]] = {"Emissivity":emissivity}

                fig, ax = plt.subplots(figsize = (12,6))
                T_min = 200
                T_max = 1000
                num_points = 100
                temperatures = np.linspace(T_min, T_max, num_points)
                plt.figure(figsize = (12,6))
                for substance_s, epsilon in table_substance.items():
                    emissivity = epsilon["Emissivity"]
                    heat_transfer_rate = emissivity * stefan_boltzmanns_constant * surface_area * temperatures ** 4
                    ax.plot(temperatures, heat_transfer_rate, label=f"{substance_s}, ε = {epsilon['Emissivity']}")
                
                #graph's info
                ax.set_xlabel("Temperature (K)")
                ax.set_ylabel("Heat Transfer Rate (W)")
                plot_title = st.text_input("Graph's title: ")
                ax.set_title(f"{plot_title}")
                ax.grid(True)
                ax.legend()
                st.pyplot(fig)

        if calculation == "8. Convective Heat Transfer":
            st.subheader("Calculation: P = hAΔT")
            A = st.number_input("Surface Area (A):", min_value = 1e-6) #unit: m²
            h = st.number_input("Convective heat transfer coefficient (h):", min_value = 0.0) #unit: W/(m²K)
            delta_T = st.number_input("Temperature difference (ΔT):", min_value = 0.0) #unit: K
            rate_of_heat_transfer = h * A *delta_T
            st.info(f"Rate of energy transfer = {rate_of_heat_transfer} W")
    
    if type_of_cal == "Greenhouse Gas Impact":
        cal_or_bar = st.selectbox("Calculate Amount of Carbon Dioxide and Methane Emitted or Compare Carbon Emissions (Bar Graph)?", ["---Please select---","Compare Carbon Emissions (Bar Graph)", "Calculate Amount of Carbon Dioxide and Methane Emitted"])
        table = [] #preparing a list for appending
        with open("list_of_fuels.csv","r", encoding = "utf-8") as file: #open csv
            reader = csv.DictReader(file)
            for row in reader:
                table.append({"Fuel": row["Fuel"],
                                "Energy Content (MJ/kg)":row["Energy Content (MJ/kg)"], "Carbon Emissions (gCO₂/MJ)":row["Carbon Emissions (gCO₂/MJ)"]})
        st.subheader("Table: Fuel, Energy Content and Carbon Emissions")
        df = pd.DataFrame(table)
        st.dataframe(df)
        st.write("*Assumes green hydrogen. No direct CO₂ emissions upon use.")
        st.write("Note: ☘️ means it is a renewable energy!")
        
        if cal_or_bar == "Calculate Amount of Carbon Dioxide and Methane Emitted":
            st.subheader("Amount of Carbon Dioxide Emitted Calculator")
            amount_of_energy = st.number_input("Based on the table above or your own value, energy content of the fuel (MJ/kg):", min_value = 0.0)
            mass = st.number_input("Mass of the fuel consumed (kg): ", min_value = 0.0)
            carbon_emissions_cal = st.number_input("Based on the table above or your own value, carbon emissions (gCO₂/MJ): ", min_value = 0.0)
            if amount_of_energy == 0 or mass == 0:
                st.warning("Energy content and Mass cannot be 0!")
                return
            carbon_emitted = amount_of_energy * mass * carbon_emissions_cal
            st.info(f"{mass}kg of the fuel emitted {round(carbon_emitted/1000,3)}kg of CO₂.")
            st.info(f"🌳 To offset {round(carbon_emitted/1000,3)}kg of CO₂, you would need approximately {round(carbon_emitted/1000/26.635)} trees!")
            methane_emissions_cal = st.number_input("Methane emissions in gCH4/MJ: ")
            methane_emitted = amount_of_energy * mass * methane_emissions_cal
            st.info(f"{mass}kg of the fuel emitted {round(methane_emitted/1000,3)}kg of CH₄.")
            st.info(f"Over a 20-year period, methane is estimated to be 80 times more potent than CO₂ when it comes to warming the planet.\n Hence, the calculated amount of methane is equivalent to {round(methane_emitted/1000 * 80, 2)}  kg of CO₂.")

        if cal_or_bar == "Compare Carbon Emissions (Bar Graph)":
            table_substance = {}
            material_names = (name["Fuel"] for name in table)
            common_substance = st.multiselect(f"Please select one or more items from the list below:", material_names) #returns list
            st.subheader("Graph")
            
            if not common_substance:
                st.warning("To display the graph, please select at least one item from the table!")
                return
            for substance in common_substance: 
                for item in table:
                    if substance == item["Fuel"]:
                        energy_content = float(item["Energy Content (MJ/kg)"])
                        carbon_emissions = float(item["Carbon Emissions (gCO₂/MJ)"])
                        table_substance[item["Fuel"]] = {"Energy Content (MJ/kg)": energy_content, "Carbon Emissions (gCO₂/MJ)": carbon_emissions}

            fig1, ax1 = plt.subplots(figsize = (12,6))
            renewable_energies = ["Hydrogen* (HHV)☘️", "Hydrogen* (LHV)☘️", "Vegetable Oil☘️", "Biodiesel☘️"] 

            fuel_names = list(table_substance.keys())
            colors = cm.get_cmap("tab20", len(fuel_names))
            fuel_color_map = {fuel: colors(i) for i, fuel in enumerate(fuel_names)}

            val_based = {k: v for k, v in sorted(table_substance.items(), key=lambda item:item[1]["Energy Content (MJ/kg)"])}
            for fuel, energy_con in val_based.items():
                Fuel = [fuel]
                Energy = [energy_con["Energy Content (MJ/kg)"]]
                energy_val = energy_con["Energy Content (MJ/kg)"]
                color = fuel_color_map[fuel]
                bars1 = ax1.bar(Fuel, Energy, width = 0.3, color=color)
                ax1.bar_label(bars1, label_type = "edge", color = "blue")
            
                if fuel in renewable_energies: #cannot use a list in ax.text!
                    ax1.text(fuel, energy_val + 5, "☘️", ha = "center", fontsize = 20)
            
            #graph's info
            ax1.set_xlabel("Fuel")
            ax1.set_ylabel("Energy Content (MJ/kg)")
            ax1.tick_params(axis = "x", labelrotation=45)
            plot_title = st.text_input("Graph's title: ", key = "F")
            ax1.set_title(f"{plot_title}")
            st.pyplot(fig1)
            
            fig2, ax2 = plt.subplots(figsize = (12,6))
            
            #table substance contains dicts within a dict, you need to specify which item in the 2nd dict to access 
            val_based = {k: v for k, v in sorted(table_substance.items(), key=lambda item:item[1]["Carbon Emissions (gCO₂/MJ)"])}
            
            for fuel, energy_con in val_based.items():
                Fuel = [fuel] #fuel gets the key, energy_con gets the value
                Energy2 = [energy_con["Carbon Emissions (gCO₂/MJ)"]]
                energy_2 = energy_con["Carbon Emissions (gCO₂/MJ)"]
                color = fuel_color_map[fuel]
                bars2 = ax2.bar(Fuel, Energy2, width = 0.3, color=color)
                ax2.bar_label(bars2, label_type = "edge", color = "blue")

                if fuel in renewable_energies: #cannot use a list in ax.text!
                    ax2.text(fuel, energy_2 + 5, "☘️", ha = "center", fontsize = 20)
            
            #graph's info
            ax2.set_xlabel("Fuel")
            ax2.set_ylabel("Carbon Emissions (gCO₂/MJ)")
            ax2.tick_params(axis = "x", labelrotation=45)
            plot_title = st.text_input("Graph's title: ", key = "C")
            ax2.set_title(f"{plot_title}")
            st.pyplot(fig2)

def main():
    type_of_calculation()

main()