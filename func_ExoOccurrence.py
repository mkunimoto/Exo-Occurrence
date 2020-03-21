import numpy as np
from scipy.stats import poisson
from scipy.stats import rayleigh
from scipy.stats import uniform
from scipy.stats import bernoulli
from scipy.stats import gamma
from scipy.stats import halfnorm

########## USER INPUT ###########

pbin = np.asarray([50.,100.])
rbin = np.asarray([1.0,1.414,2.0,2.828,4.0,5.657])
stars = np.loadtxt("FGK_properties.dat")

########### CONSTANTS ###########

earth_radius_in_m = 6.378e6
sun_radius_in_m = 6.9566e8
rsol_in_au = 0.00465
rearth_in_rsol = 0.009168

kepler_exp_time_internal = 6.019802903/(24.*60.*60.)
kepler_read_time_internal = 0.5189485261/(24.*60.*60.)
num_exposures_per_LC = 270.
LC_integration_time = kepler_exp_time_internal*num_exposures_per_LC
LC_read_time = kepler_read_time_internal*num_exposures_per_LC
LC_duration = LC_integration_time + LC_read_time
LC_rate = 1./LC_duration

Nt = len(stars)

########### FUNCTIONS ###########

def semimajor_axis_radius_ratio(P, M, R):
	a3 = M*(P/365.)**2.
	a = a3**(1./3.)
	return a/R/rsol_in_au

def transit_depth(Rp, R, u1, u2):
	k = Rp/R*rearth_in_rsol
	c0 = 1. - (u1 + u2)
	w = c0/4. + (u1 + 2.*u2)/6. - u2/8.
	return 1. - (c0/4. + (u1+2.*u2)*(1.-k**2.)**(3./2.)/6. - u2*(1.-k**2)/8.)/w

def transit_duration(P, aR, b):
        return 24.*(P/np.pi)/aR*np.sqrt(1.-b**2.)

def impact_parameter(aR, i):
	return aR*np.cos(i)

def num_transits(P, Tobs, fduty):
	return np.round(Tobs/P*fduty)

def transit_SNR(P, Rp, dur, R, u1, u2, Tobs, fduty, std):
	Ntr = num_transits(P, Tobs, fduty)
	depth = transit_depth(Rp, R, u1, u2)
	return depth/std*np.sqrt(Ntr*dur*LC_rate)

def P_win(P, Tobs, fduty):
	j = Tobs/P
	return np.where(j < 2., 0.0, 1.-(1.-fduty)**j-j*fduty*(1.-fduty)**(j-1.)-0.5*j*(j-1.)*fduty*fduty*(1.-fduty)**(j-2.))

def det_comp(SNR, A, B, C):
        return C*gamma.cdf(SNR, A, scale=B)

def P_det(Ntr, SNR):
        return np.where(Ntr >= 37, det_comp(SNR,12.23315232,  0.78203581,  0.94645662), np.where(Ntr >= 19, det_comp(SNR,14.86511928,  0.72917663,  0.91642721), np.where(Ntr >= 10, det_comp(SNR,11.45382365,  1.07249528,  0.91176524), np.where(Ntr >= 7, det_comp(SNR,11.54128644,  1.20628098,  0.88169029), np.where(Ntr == 6, det_comp(SNR,11.48116478,  1.27632116,  0.83694848), np.where(Ntr == 5, det_comp(SNR,13.54878807,  1.09003313,  0.7704247), np.where(Ntr == 4, det_comp(SNR,17.47440559,  0.90589395,  0.66508744), np.where(Ntr == 3, det_comp(SNR,12.02912833,  1.38916308,  0.46525859), 0.0))))))))

def transit_noise_model(P, Rp, R, aR, b, Tobs, fduty, std):
	t0 = P*np.random.uniform()
	tao0 = P/(aR*2.*np.pi)
	r = rearth_in_rsol*Rp/R
	T = 2.*tao0*np.sqrt(1.-b**2.)
	tau = 2.*tao0*r/np.sqrt(1.-b**2.)
	Ttot = P
	I = LC_integration_time
	Lambda = LC_rate*num_transits(P, Tobs, fduty)
	tau3 = tau**3.
	I3 = I**3.
	a2 = (5.*tau3 + I3 - 5.*tau*tau*I)/tau3
	a3 = (9.*I**5.*Ttot - 40.*tau3*I*I*Ttot + 120.*tau**4.*I*(3.*Ttot - 2.*tau))/tau**6.
	a4 = (a3*tau**5. + I**4.*(54.*tau - 35.*Ttot) - 12.*tau*I3*(4.*tau+Ttot) + 360.*tau**4.*(tau - Ttot))/tau**5.
	a5 = (a2*(24.*T*T*(I-3.*tau) - 24.*T*Ttot*(I-3.*tau)) + tau3*a4)/tau3
	a11 = (I*Ttot - 3.*tau*(Ttot - 2.*tau))/(tau*tau)
	b1 = (6.*I*I - 3.*I*Ttot + tau*Ttot)/(I*I)
	b4 = (6.*T*T - 6.*T*Ttot + I*(5.*Ttot - 4.*I))/(I*I)
	b6 = (12.*b4*I3 + 4.*tau*(-6.*T*T + 6.*T*Ttot + I*(13.*Ttot - 30.*I)))/I3
	b7 = (b6*I**5. + 4.*tau*tau*I*I*(12.*I - 11.*Ttot) + tau3*I*(11.*Ttot - 6.*I) - tau**4.*Ttot)/I**5.
	depth = r**2.
	sigma_depth = np.where(tau >= I, np.sqrt(abs(24.*a11*a2/(tau*a5))), np.sqrt(abs(24.*b1/(I*b7))))
	sigma_depth *= std/np.sqrt(Lambda)
	return np.random.normal(loc=depth, scale=sigma_depth)

def assumed_stellar_radius(sigdown, sigup, flag):
	return np.where(flag == 1, halfnorm.rvs(scale=sigup), -halfnorm.rvs(scale=sigdown))

####################

def simulate_catalogue(params):
        rates = [params['r1'],params['r2'],params['r3'],params['r4'],params['r5']] 
	total_rate = np.sum(rates)
	prob_bin = rates/total_rate
	num_bin = len(rates)
	num_pl = poisson.rvs(total_rate*Nt)
	if (num_bin == 1):
		allRp = np.exp(np.random.uniform(np.log(rbin[0]),np.log(rbin[1]),size=num_pl))
	else:
		bins = np.random.choice(num_bin, num_pl, p=prob_bin)
		allRp = np.exp(np.random.uniform(np.log(rbin[bins]),np.log(rbin[bins+1])))
	allP = np.exp(np.random.uniform(np.log(pbin[0]),np.log(pbin[1]),size=num_pl))
	allID = np.random.random_integers(0,high=len(stars)-1,size=num_pl)
        # Generate orbital properties
	alli = np.arccos(np.random.uniform(size=num_pl))
	allR = stars[allID,1]
	allaR = semimajor_axis_radius_ratio(allP, stars[allID,5], allR)
	allb = impact_parameter(allaR, alli)
        # Remove non-transiting planets
	itrans = (allb <= 1.0)
	transP = allP[itrans]
	transRp = allRp[itrans]
	transR = allR[itrans]
        # Calculate SNR and detection probabilities
	dur = transit_duration(transP, allaR[itrans], allb[itrans])
        Ntr = num_transits(transP, stars[allID,12][itrans], stars[allID,13][itrans])
	SNR = transit_SNR(transP, transRp, dur, transR, stars[allID,9][itrans], stars[allID,10][itrans], stars[allID,12][itrans], stars[allID,13][itrans], stars[allID,11][itrans])
	Pdet = P_det(Ntr, SNR)
	Pwin = P_win(transP, stars[allID,12][itrans], stars[allID,13][itrans])
        # Remove non-detected planets
	idet = (bernoulli.rvs(Pdet*Pwin) == 1)
	detRp = transRp[idet]
        detP = transP[idet]
        detR = transR[idet]
	# Generated observed properties
	obsdepth = transit_noise_model(detP,detRp,detR,allaR[itrans][idet],allb[itrans][idet],stars[allID,12][itrans][idet],stars[allID,13][itrans][idet],stars[allID,11][itrans][idet])
	obsdepth[obsdepth < 0.0] = 0.0
	flag = np.random.randint(0,2,size=len(detRp))
	obsR = detR + assumed_stellar_radius(stars[allID,2][itrans][idet],stars[allID,3][itrans][idet],flag)
	obsRp = np.sqrt(obsdepth)*obsR/rearth_in_rsol
	data_sim = np.column_stack((detP, obsRp))
	return data_sim

def distance(data_obs, params):
	data_sim = params['dataset1']
        per_sim = data_sim[:,0]
        per_obs = data_obs[:,0]
        Rp_sim = data_sim[:,1]
	Rp_obs = data_obs[:,1]
	num_bin = len(rbin) - 1
	sum_obs = np.zeros(num_bin)
	sum_sim = np.zeros(num_bin)
	for i in range(num_bin):
		sum_obs[i] = len(Rp_obs[(per_obs >= pbin[0]) & (per_obs <= pbin[1]) & (Rp_obs >= rbin[i]) & (Rp_obs <= rbin[i+1])])
		sum_sim[i] = len(Rp_sim[(per_sim >= pbin[0]) & (per_sim <= pbin[1]) & (Rp_sim >= rbin[i]) & (Rp_sim <= rbin[i+1])])
	sum_obs = sum_obs/Nt
	sum_sim = sum_sim/Nt
	dist = []
	for i in range(num_bin):
		if (sum_obs[i] == sum_sim[i]):
			dist.append(0.0)
		else:
			dist.append(np.abs(sum_obs[i]-sum_sim[i])/np.sqrt(sum_obs[i]+sum_sim[i]))
	return np.atleast_1d(np.sum(dist))

