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
                  SeparatedFlow, SuctionSide, PressureSide, Aoa, ChordBasedReynolds, Transition, rot_speed, u_wind_speed, TSRtd, dstar_type, state_ss_us, anglesteps, phi_type, theta_type, IntegralLengthScale, TurbulenceIntensity,rot_speed_check, u_wind_speed_check, TSR_check, Lowson_type, obs_x_pos, obs_y_pos, obs_z_pos, obs_x_pos_rotor, obs_y_pos_rotor, obs_z_pos_rotor, tower_height, initial_azimuth, yaw_angle, tower_to_hub_distance, tower_to_rotor_distance, rotation_type, number_loops, time, timesteps, shear_roughness, shear_height, shear_speed, shear_check, qs3DSim, valRel_TE_check, valRel_TE, valReu_TE_check, valReu_TE, valMal_TE_check, valMal_TE, valMau_TE_check, valMau_TE, valAOAl_TE_check, valAOAl_TE, valAOAu_TE_check, valAOAu_TE, valRel_LE_check, valRel_LE, valReu_LE_check, valReu_LE, valMal_LE_check, valMal_LE, valMau_LE_check, valMau_LE}; //Sara and Alexandre MOD
	};
};

#endif // PARAMETERKEYS_H
