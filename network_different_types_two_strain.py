                    # if (S[contact] == 1):
                    #     if random.random() < beta * dt:
                    #         I2[contact] = 1
                    #         S[contact] = 0
                    # elif (R1[contact] == 1 and V[contact] == 0) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # Recovered/protected; if just recovered from strain 1, transmission probability is NOT reduced.
                    #     effective_beta = beta 
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1
                    # elif (V[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # If vaccinated, transmission probability IS reduced.
                    #     effective_beta = beta * (1.0 - vaccine_effectiveness)
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1

                    # ## 2. Vaccine not universal, recovered different for the two strains
                    # if (S[contact] == 1):
                    #     if random.random() < beta * dt:
                    #         I2[contact] = 1
                    #         S[contact] = 0
                    # elif (R1[contact] == 1 or V[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # Recovered/protected; if just recovered from strain 1, transmission probability is NOT reduced.
                    #     effective_beta = beta
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1

                    # ## 3. Vaccine universal, recovered same for the two strains
                    # if (S[contact] == 1):
                    #     if random.random() < beta * dt:
                    #         I2[contact] = 1
                    #         S[contact] = 0

                    # elif (V[contact] == 1 or R1[contact]) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # If vaccinated, transmission probability IS reduced.
                    #     effective_beta = beta * (1.0 - vaccine_effectiveness)
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1

                    ## 4. Vaccine not universal, recovered same for the two strains
                    if (S[contact] == 1):
                        if random.random() < beta * dt:
                            I2[contact] = 1
                            S[contact] = 0
                    elif (V[contact] == 1 and R1[contact] == 0) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                        # Vaccinated; if just vaccinated , transmission probability is NOT reduced.
                        effective_beta = beta 
                        if random.random() < effective_beta * dt:
                            R1[contact] = 0
                            I2[contact] = 1
                    elif (R1[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                        # If recpvered, transmission probability IS reduced.
                        effective_beta = beta * (1.0 - vaccine_effectiveness)
                        if random.random() < effective_beta * dt:
                            R1[contact] = 0
                            I2[contact] = 1