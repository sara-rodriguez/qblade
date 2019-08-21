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
		enum Key {Name, Rho, Viscosity, Discretize, MaxIterations, MaxEpsilon, RelaxFactor, TipLoss, VariableInduction,
				  TipSpeedFrom, TipSpeedTo, TipSpeedDelta, Windspeed};
	};

	struct CDMSData {
		enum Key {Name, Rho, Viscosity, Discretize, MaxIterations, MaxEpsilon, RelaxFactor, TipLoss, VariableInduction,
				  WindspeedFrom, WindspeedTo, WindspeedDelta, RotationalFrom, RotationalTo, RotationalDelta, PitchFrom,
				  PitchTo, PitchDelta};
	};
	struct TDMSData {
		enum Key {Name, Rho, Viscosity, Discretize, MaxIterations, MaxEpsilon, RelaxFactor, WindspeedFrom, WindspeedTo,
				  WindspeedDelta, TipLoss, VariableInduction, AnnualYield};
	};
	struct NoiseSimulation {
        enum Key {Name, WettedLength, DistanceObsever, OriginalVelocity, OriginalChordLength, OriginalMach,
                  DStarChordStation, DStarScalingFactor, EddyConvectionMach, DirectivityTheta, DirectivityPhi,
                  SeparatedFlow, SuctionSide, PressureSide, Aoa, ChordBasedReynolds, Transition, sects, rot_speed, u_wind_speed, TSRtd, dstar_type, phi_type, IntegralLengthScale, TurbulenceIntensity,rot_speed_calc, u_wind_speed_calc, TSR_calc,rot_speed_check, u_wind_speed_check, TSR_check}; //Sara and Alexandre MOD
	};
};

#endif // PARAMETERKEYS_H
