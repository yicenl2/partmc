CDF       
      gas_species    X   aero_species      aero_source       aero_weight_group         aero_weight_class         aero_particle         aero_removed            title          PartMC version 2.6.1 output file   source        PartMC version 2.6.1   UUID      $39618B07-1B6F-4A9A-B351-8623BF88457B   history       =2023-06-14T09:27:38.530-05:00 created by PartMC version 2.6.1      Conventions       CF-1.4        /   time             unit      s      description       #time elapsed since simulation start          l   timestep             unit      s      description       current timestep size            t   timestep_index               description       Jan integer that is 1 on the first timestep, 2 on the second timestep, etc.           |   repeat               description       2repeat number of this simulation (starting from 1)           �   temperature              unit      K      standard_name         air_temperature          �   relative_humidity                unit      1      standard_name         relative_humidity            �   pressure             unit      Pa     standard_name         air_pressure         �   	longitude                unit      degree_east    standard_name         	longitude            �   latitude             unit      degree_north   standard_name         latitude         �   altitude             unit      m      standard_name         altitude         �   start_time_of_day                unit      s      description       9time-of-day of simulation start in seconds since midnight            �   start_day_of_year                description       &day-of-year number of simulation start           �   elapsed_time             unit      s      description       #elapsed time since simulation start          �   solar_zenith_angle               unit      radian     description       (current angle from the zenith to the sun         �   height               unit      m      	long_name         boundary layer mixing height         �   gas_species                 names        �NO2,NO,O,O3,NO3,O1D,OH,HO2,N2O5,HNO3,HONO,PNA,H2O2,XO2,XO2N,NTR,ROOH,FORM,ALD2,ALDX,PAR,CO,MEO2,MEPX,MEOH,HCO3,FACD,C2O3,PAN,PACD,AACD,CXO3,PANX,ROR,OLE,ETH,IOLE,TOL,CRES,TO2,TOLRO2,OPEN,CRO,MGLY,XYL,XYLRO2,ISOP,ISPD,ISOPRXN,TERP,TRPRXN,SO2,SULF,SULRXN,ETOH,ETHA,CL2,CL,HOCL,CLO,FMCL,HCL,TOLNRXN,TOLHRXN,XYLNRXN,XYLHRXN,BENZENE,BENZRO2,BNZNRXN,BNZHRXN,SESQ,SESQRXN,M,O2,N2,H2O,CH4,H2,N2O,DUMMY,NH3,DMS,ISOP-P1,ISOP-P2,TERP-P1,TERP-P2,COV,APIN     description       tdummy dimension variable (no useful value) - read species names as comma-separated values from the 'names' attribute     `   �   gas_mosaic_index                	long_name         MOSAIC indices of gas species        `  "8   gas_mixing_ratio                unit      ppb    	long_name         mixing ratios of gas species     �  #�   aero_species               names        �dust.LHD_DUST,dust.LLD_DUST,sea_salt.SEA_SALT,organic_matter.POA,organic_matter.SOA,organic_matter.COV_aero,organic_matter.COV-PROD_aero,black_carbon.BC_phob,black_carbon.BC_phil,inorganics_dry.AS,inorganics_dry.COV_aero,inorganics_dry.COV-PROD_aero,inorganics_wet.AS,inorganics_wet.H2O_aq,inorganics_wet.COV_aero,inorganics_wet.COV-PROD_aero,other_PM.other_PM,other_PM.other_other_PM,other_PM.COV_aero,other_PM.COV-PROD_aero      description       tdummy dimension variable (no useful value) - read species names as comma-separated values from the 'names' attribute      P  &X   aero_source                names         dry_ammonium_sulfate   description       sdummy dimension variable (no useful value) - read source names as comma-separated values from the 'names' attribute         &�   aero_mosaic_index                  	long_name         !MOSAIC indices of aerosol species         P  &�   aero_density               unit      kg/m^3     	long_name         densities of aerosol species      �  &�   aero_num_ions                  	long_name         4number of ions after dissociation of aerosol species      P  '�   aero_molec_weight                  unit      kg/mol     	long_name         $molecular weights of aerosol species      �  '�   
aero_kappa                 unit      1      	long_name         5hygroscopicity parameters (kappas) of aerosol species         �  (�   fractal_dimension                unit      1      description       !particle volume fractal dimension           ),   fractal_prime_radius             unit      m      description       radius of primary particles         )4   fractal_vol_fill_factor              unit      1      description       volume filling factor           )<   aero_weight_group                  description       *dummy dimension variable (no useful value)          )D   aero_weight_class                  description       *dummy dimension variable (no useful value)          )H   weight_type                   description       �type of each aerosol weighting function: 0 = invalid, 1 = none (w(D) = 1), 2 = power (w(D) = (D/D_0)^alpha), 3 = MFA (mass flow) (w(D) = (D/D_0)^(-3))          )L   weight_magnitude                  unit      m^{-3}     description       %magnitude for each weighting function           )P   weight_exponent                   unit      1      description       Oexponent alpha for the power weight_type, set to -3 for MFA, and zero otherwise         )X   aero_particle                  description       *dummy dimension variable (no useful value)          )`   aero_particle_mass                    unit      kg     	long_name         +constituent masses of each aerosol particle      �  )x   aero_n_orig_part                  	long_name         gnumber of original constituent particles from each source that coagulated to form each aerosol particle         -8   aero_particle_weight_group                 	long_name         ,weight group number of each aerosol particle        -P   aero_particle_weight_class                 	long_name         ,weight class number of each aerosol particle        -h   aero_water_hyst_leg                	long_name         >leg of the water hysteresis curve leg of each aerosol particle          -�   aero_num_conc                  unit      m^{-3}     	long_name         &number concentration for each particle        0  -�   aero_id                	long_name         )unique ID number of each aerosol particle           -�   aero_least_create_time                 unit      s      	long_name         ,least creation time of each aerosol particle   description       �least (earliest) creation time of any original constituent particles that coagulated to form each particle, measured from the start of the simulation         0  -�   aero_greatest_create_time                  unit      s      	long_name         /greatest creation time of each aerosol particle    description       �greatest (latest) creation time of any original constituent particles that coagulated to form each particle, measured from the start of the simulation        0  .   aero_removed               description       *dummy dimension variable (no useful value)          .@   aero_removed_id                	long_name         ID of removed particles         .D   aero_removed_action                	long_name         reason for particle removal    description       �valid is 0 (invalid entry), 1 (removed due to dilution), 2 (removed due to coagulation -- combined particle ID is in \c aero_removed_other_id), 3 (removed due to populating halving), or 4 (removed due to weighting changes           .H   aero_removed_other_id                  	long_name         (ID of other particle involved in removal   description       �if <tt>aero_removed_action(i)</tt> is 2 (due to coagulation), then <tt>aero_removed_other_id(i)</tt> is the ID of the resulting combined particle, or 0 if the new particle was not created         .L@�`     @$         5   @r�fffff        @�j             @V�             @�        �@�`             ?�                                 	   
                                                                      !   "   #   $   %   &   '   (   )   *   +   ,   -   .   /   0   1   2   3   4   5   6   7   8   9   :   ;   <   =   >   ?   @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X                                                                                                                                                                                                                                                                                                                                                                @$7�=�O@ ����J>����{�@HB_��4�?0UZ(OI=E|Q���? �5 ?piS�Q��?Uۘ��;@?�փ��xj?�K�����?�sm��\�?�xn@HP?b��o��*?!W3�_�?�ΞK٪?�y&u(@u�l�v?��{5G��?ڥ+0w	@A�2��8@v���Ğo?62�ѐ{f?ޗ�.U�?�BzQ�>�����?s˃J@�?:���e6?�!WՂ?X�"�.Q?ɍ��e>�ňF#��?�f'���>HY�Yk��?��=�@I�+��?��O2�?�!z�?$S?��̻a�k>� >���p>���^�n�?c8�S�>��<Q/r�?��N�a?�hr��#�?� F�@��PN�?�i�� z@8�H��A4
�e|ґ4&x���@T�\0h{?ͦ���?ͦ���4Y�!`.?���O�B=;`3O%C4>Deg7c>���>�.|>���
ў�?[A���u?�X� ��?�����g?&xzv�h_?�P���?F�`%�84���c3�m����N6A�܏Q�45�$�����4��d�Y�4o����.A��e    A��m�2�(A�>�   4�����@�/Ւq
@�҃ ��4����A?��mXQ�@��lG�j�4����A4����A4����A4����A4����A=��`�4����A                           	   
                                 "y  +�"y  +�                ��0    SUD_   T           !    ��@    ��P    @��     @��     @��     @�      @�      @��     @��     @�@     @�@     @��     @��     @��     @��     @�@     @��     @��     @�      @�      @��     @��                                                                                   ?�������?�������?�������?ڋC��%?����^?ə�����?�(�\)?�������?�������?����?)?ə�����?�(�\)?����?)?�r�����?ə�����?�(�\)?�������?�������?ə�����?�(�\)?�������?�������?��\(�        ?�������?�������?�������        ?PbM���?�������?�������?�������?�������        ?�������?�������                ?�������?�������@      >Ey��0�:?�               A�ׄ                              3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����%3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&<s`���<s`���<s`���<s`���<s`���<s`���9���+��9���+��9���+��9���+��9���+��9���+��3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&;�{�[�;�{�[�;�{�[�;�{�[�;�{�[�;�{�[�<s`���H<s`���H<s`���H<s`���H<s`���H<s`���H9�n(��j9�n(��j9�n(��j9�n(��j9�n(��j9�n(��j3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&<s`���<s`���<s`���<s`���<s`���<s`���9�9�٘)9�9�٘)9�9�٘)9�9�٘)9�9�٘)9�9�٘)3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&3y����&<s`���<s`���<s`���<s`���<s`���<s`���9�X���9�X���9�X���9�X���9�X���9�X���                                                                              A�ׄ    A�ׄ    A�ׄ    A�ׄ    A�ׄ    A�ׄ                                                                                                                                     