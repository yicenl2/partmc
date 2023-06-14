CDF       
      gas_species    X   aero_species      aero_source       aero_weight_group         aero_weight_class         aero_particle         aero_removed            title          PartMC version 2.6.1 output file   source        PartMC version 2.6.1   UUID      $8F8F78A5-F9CE-4744-B747-8CE7BA3FD16A   history       =2023-06-14T13:55:28.427-05:00 created by PartMC version 2.6.1      Conventions       CF-1.4        /   time             unit      s      description       #time elapsed since simulation start          @   timestep             unit      s      description       current timestep size            H   timestep_index               description       Jan integer that is 1 on the first timestep, 2 on the second timestep, etc.           P   repeat               description       2repeat number of this simulation (starting from 1)           T   temperature              unit      K      standard_name         air_temperature          X   relative_humidity                unit      1      standard_name         relative_humidity            `   pressure             unit      Pa     standard_name         air_pressure         h   	longitude                unit      degree_east    standard_name         	longitude            p   latitude             unit      degree_north   standard_name         latitude         x   altitude             unit      m      standard_name         altitude         �   start_time_of_day                unit      s      description       9time-of-day of simulation start in seconds since midnight            �   start_day_of_year                description       &day-of-year number of simulation start           �   elapsed_time             unit      s      description       #elapsed time since simulation start          �   solar_zenith_angle               unit      radian     description       (current angle from the zenith to the sun         �   height               unit      m      	long_name         boundary layer mixing height         �   gas_species                 names        �NO2,NO,O,O3,NO3,O1D,OH,HO2,N2O5,HNO3,HONO,PNA,H2O2,XO2,XO2N,NTR,ROOH,FORM,ALD2,ALDX,PAR,CO,MEO2,MEPX,MEOH,HCO3,FACD,C2O3,PAN,PACD,AACD,CXO3,PANX,ROR,OLE,ETH,IOLE,TOL,CRES,TO2,TOLRO2,OPEN,CRO,MGLY,XYL,XYLRO2,ISOP,ISPD,ISOPRXN,TERP,TRPRXN,SO2,SULF,SULRXN,ETOH,ETHA,CL2,CL,HOCL,CLO,FMCL,HCL,TOLNRXN,TOLHRXN,XYLNRXN,XYLHRXN,BENZENE,BENZRO2,BNZNRXN,BNZHRXN,SESQ,SESQRXN,M,O2,N2,H2O,CH4,H2,N2O,DUMMY,NH3,DMS,ISOP-P1,ISOP-P2,TERP-P1,TERP-P2,COV,APIN     description       tdummy dimension variable (no useful value) - read species names as comma-separated values from the 'names' attribute     `   �   gas_mosaic_index                	long_name         MOSAIC indices of gas species        `  "   gas_mixing_ratio                unit      ppb    	long_name         mixing ratios of gas species     �  #l   aero_species               names        �dust.LHD_DUST,dust.LLD_DUST,sea_salt.SEA_SALT,organic_matter.POA,organic_matter.SOA,organic_matter.COV_aero,organic_matter.COV-PROD_aero,black_carbon.BC_phob,black_carbon.BC_phil,inorganics_dry.AS,inorganics_dry.COV_aero,inorganics_dry.COV-PROD_aero,inorganics_wet.AS,inorganics_wet.H2O_aq,inorganics_wet.COV_aero,inorganics_wet.COV-PROD_aero,other_PM.other_PM,other_PM.other_other_PM   description       tdummy dimension variable (no useful value) - read species names as comma-separated values from the 'names' attribute      H  &,   aero_source                names         dry_ammonium_sulfate   description       sdummy dimension variable (no useful value) - read source names as comma-separated values from the 'names' attribute         &t   aero_mosaic_index                  	long_name         !MOSAIC indices of aerosol species         H  &x   aero_density               unit      kg/m^3     	long_name         densities of aerosol species      �  &�   aero_num_ions                  	long_name         4number of ions after dissociation of aerosol species      H  'P   aero_molec_weight                  unit      kg/mol     	long_name         $molecular weights of aerosol species      �  '�   
aero_kappa                 unit      1      	long_name         5hygroscopicity parameters (kappas) of aerosol species         �  ((   fractal_dimension                unit      1      description       !particle volume fractal dimension           (�   fractal_prime_radius             unit      m      description       radius of primary particles         (�   fractal_vol_fill_factor              unit      1      description       volume filling factor           (�   aero_weight_group                  description       *dummy dimension variable (no useful value)          (�   aero_weight_class                  description       *dummy dimension variable (no useful value)          (�   weight_type                   description       �type of each aerosol weighting function: 0 = invalid, 1 = none (w(D) = 1), 2 = power (w(D) = (D/D_0)^alpha), 3 = MFA (mass flow) (w(D) = (D/D_0)^(-3))          (�   weight_magnitude                  unit      m^{-3}     description       %magnitude for each weighting function           (�   weight_exponent                   unit      1      description       Oexponent alpha for the power weight_type, set to -3 for MFA, and zero otherwise         (�   aero_particle                  description       *dummy dimension variable (no useful value)           (�   aero_particle_mass                    unit      kg     	long_name         +constituent masses of each aerosol particle      �  )   aero_n_orig_part                  	long_name         gnumber of original constituent particles from each source that coagulated to form each aerosol particle          -�   aero_particle_weight_group                 	long_name         ,weight group number of each aerosol particle         -�   aero_particle_weight_class                 	long_name         ,weight class number of each aerosol particle         -�   aero_water_hyst_leg                	long_name         >leg of the water hysteresis curve leg of each aerosol particle           -�   aero_num_conc                  unit      m^{-3}     	long_name         &number concentration for each particle        @  .   aero_id                	long_name         )unique ID number of each aerosol particle            .L   aero_least_create_time                 unit      s      	long_name         ,least creation time of each aerosol particle   description       �least (earliest) creation time of any original constituent particles that coagulated to form each particle, measured from the start of the simulation         @  .l   aero_greatest_create_time                  unit      s      	long_name         /greatest creation time of each aerosol particle    description       �greatest (latest) creation time of any original constituent particles that coagulated to form each particle, measured from the start of the simulation        @  .�   aero_removed               description       *dummy dimension variable (no useful value)          .�   aero_removed_id                	long_name         ID of removed particles         .�   aero_removed_action                	long_name         reason for particle removal    description       �valid is 0 (invalid entry), 1 (removed due to dilution), 2 (removed due to coagulation -- combined particle ID is in \c aero_removed_other_id), 3 (removed due to populating halving), or 4 (removed due to weighting changes           .�   aero_removed_other_id                  	long_name         (ID of other particle involved in removal   description       �if <tt>aero_removed_action(i)</tt> is 2 (due to coagulation), then <tt>aero_removed_other_id(i)</tt> is the ID of the resulting combined particle, or 0 if the new particle was not created         .�@�      @$         !   @r�fffff        @�j             @V�             @�        �@�              ?�                                 	   
                                                                      !   "   #   $   %   &   '   (   )   *   +   ,   -   .   /   0   1   2   3   4   5   6   7   8   9   :   ;   <   =   >   ?   @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X                                                                                                                                                                                                                                                                                                                                                                @���R�?��;�>��� 	�@Hi��4�?.��!I�=E�%C0p ?zVD~�o?wg��c�?J���{/=?�k�6Z�2?��_Gˍ?��-N*�?��i�?mꌮ��?/�+����?Ї�WƵ\?�;���=�@A<���?�?�n��?�|I��@6V��U@r����K?@(���:�?�2WD|Q?�醽��>��"��,?`1��t�?!+����$?��g��g?VI �H_�?ɘ�. �>�_w�?��Pi�>@�+h'~?�*� �@�o��=�?�@F��?�,$��g�?��2A�&�>���Y;u�>����ߑ�?J<�q�v>��0�C�?��bn{B\?�@pV�?^�NS@�8\�?������^?���x4
����`4"vp���@I��x�#k?��d�9?��d�94S����?���>yO}=:*-~�>JE
�z>���床4>�ܢ/L ?Q�f��}?�]���K?�<�C7�?�cğ[�?��(.9�?@4�3c4�z:�3��ᬉ��4��B �40�� �48��;x4bbͷpqA��e    A��nI}�A�>�   4�����@�/����+@���734����A?ᢥ��@z��D� �4����A4����A4����A4����A4����A@H��2Ȉ4����A                           	   
                           ���  *����  *�                �H     m_ciettaOC.rea_V  or            @��     @��     @��     @�      @�      @��     @��     @�@     @�@     @��     @��     @��     @��     @�@     @��     @��     @�      @�                                                                            ?�������?�������?�������?ڋC��%?����^?ə�����?ə�����?�������?�������?����?)?ə�����?ə�����?����?)?�r�����?ə�����?ə�����?�������?�������?�������?�������?��\(�        ?�������?�������?�������        ?PbM���?�������?�������?�������?�������        ?�������?�������                @      >Ey��0�:?�               A�ׄ                                    3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&9I���˽�9I���˽�9I���˽�9I���˽�9I���˽�9I���˽�9I���˽�9I���˽�3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&;�{�[�;�{�[�;�{�[�;�{�[�;�{�[�;�{�[�;�{�[�;�{�[�:u�>{��`:u�>{��`:u�>{��`:u�>{��`:u�>{��`:u�>{��`:u�>{��`:u�>{��`3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&;�'a���;�'a���;�'a���;�'a���;�'a���;�'a���;�'a���;�'a���3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&                                                                                                        A�ׄ    A�ׄ    A�ׄ    A�ׄ    A�ׄ    A�ׄ    A�ׄ    A�ׄ                                                                                                                                                                           