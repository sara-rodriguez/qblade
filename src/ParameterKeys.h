#ifndef PARAMETERKEYS_H
#define PARAMETERKEYS_H


/* This enum approach is preferred over a string approach, because the misspelling of a key will be cought by the
 * compiler. Additionally this approach should perform better, though it is insignificant.
 * */

struct Parameter {
	struct Windfield {
		enum Key {Name, Time, Timesteps, Points, RotorRadius, HubHeight, WindSpeed, Turbulence, ShearLayer,
				  MeasurementHeight, RoughnessLength};
	};
	struct TData {
		enum Key {Name, Blade, VCutIn, VCutOut, TurbineOffset, TurbineHeight, RotorHeight, RotorMaxRadius,
				  RotorSweptArea, LossFactor, FixedLosses, RotationalSpeed, RotationalSpeedMin, RotationalSpeedMax, TSR};
	};
	struct DMSData {
        enum Key {Name, Temperature, Rho, Viscosity, Discretize, MaxIterations, MaxEpsilon, RelaxFactor, TipLoss, VariableInduction,
                  TipSpeedFrom, TipSpeedTo, TipSpeedDelta, Windspeed};//Sara
	};

	struct CDMSData {
        enum Key {Name, Temperature, Rho, Viscosity, Discretize, MaxIterations, MaxEpsilon, RelaxFactor, TipLoss, VariableInduction,
				  WindspeedFrom, WindspeedTo, WindspeedDelta, RotationalFrom, RotationalTo, RotationalDelta, PitchFrom,
                  PitchTo, PitchDelta};//Sara
	};
	struct TDMSData {
        enum Key {Name, Temperature, Rho, Viscosity, Discretize, MaxIterations, MaxEpsilon, RelaxFactor, WindspeedFrom, WindspeedTo,
                  WindspeedDelta, TipLoss, VariableInduction, AnnualYield};//Sara
	};
	struct NoiseSimulation {
        enum Key {Name, WettedLength, DistanceObsever, OriginalVelocity, OriginalChordLength, OriginalMach,
                  DStarChordStation, DStarScalingFactor, EddyConvectionMach, DirectivityTheta, DirectivityPhi,
                  SeparatedFlow, SuctionSide, PressureSide, Aoa, ChordBasedReynolds, Transition, rot_speed, u_wind_speed, TSRtd, dstar_type, state_ss_us, anglesteps, phi_type, theta_type, IntegralLengthScale, TurbulenceIntensity,rot_speed_check, u_wind_speed_check, TSR_check, Lowson_type, obs_x_pos, obs_y_pos, obs_z_pos, obs_x_pos_rotor, obs_y_pos_rotor, obs_z_pos_rotor, tower_height, initial_azimuth, yaw_angle, tower_to_hub_distance, rotation_type, number_loops, time, timesteps, shear_roughness, shear_height, shear_speed, shear_check, qs3DSim, valRel_TE_check, valRel_TE, valReu_TE_check, valReu_TE, valMal_TE_check, valMal_TE, valMau_TE_check, valMau_TE, valAOAl_TE_check, valAOAl_TE, valAOAu_TE_check, valAOAu_TE, valRel_LBL_VS_check, valRel_LBL_VS, valReu_LBL_VS_check, valReu_LBL_VS, valMal_LBL_VS_check, valMal_LBL_VS, valMau_LBL_VS_check, valMau_LBL_VS, valAOAl_LBL_VS_check, valAOAl_LBL_VS, valAOAu_LBL_VS_check, valAOAu_LBL_VS, valRel_LE_check, valRel_LE, valReu_LE_check, valReu_LE, valMal_LE_check, valMal_LE, valMau_LE_check, valMau_LE, autopolars_check, LBLVS, valRel_tipvortex_check, valRel_tipvortex, valReu_tipvortex_check, valReu_tipvortex, valMal_tipvortex_check, valMal_tipvortex, valMau_tipvortex_check, valMau_tipvortex, valAOAl_tipvortex_check, valAOAl_tipvortex, valAOAu_tipvortex_check, valAOAu_tipvortex, blunt_check, tipvortex_check,  valRel_blunt_check, valRel_blunt, valReu_blunt_check, valReu_blunt, valMal_blunt_check, valMal_blunt, valMau_blunt_check, valMau_blunt, valAOAl_blunt_check, valAOAl_blunt, valAOAu_blunt_check, valAOAu_blunt, hblunt_check, hblunt,flat_tip_check, valPsil_check, valPsiu_check, valPsil, valPsiu, vegetation, propagation_check, rel_humidity, vegetation_check, atm_check, LE_check, qs3d_check}; //Sara and Alexandre MOD
	};
};

#endif // PARAMETERKEYS_H
