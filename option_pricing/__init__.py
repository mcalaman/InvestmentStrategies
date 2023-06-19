import math
from scipy.stats import norm
import pandas as pd
import os
import datetime
import yfinance as yf


class BSOption:
    """ General cloass to price BS options"""
    def __init__(self, S, K, r, q, T, sigma):
        self.S = S
        self.K = K
        self.r = r
        self.q = q
        self.T = T
        self.sigma = sigma
        
        self.d1 = (math.log(self.S / self.K) + (self.r - self.q + (self.sigma ** 2) / 2) * self.T) / (self.sigma * math.sqrt(self.T))
        self.d2 = self.d1 - self.sigma * math.sqrt(self.T)
        #print(self.__class__.__name__)


class BSEuropeanPut(BSOption):
    """ BS Put options"""
    def __init__(self, S, K, r, q, T, sigma):
        BSOption.__init__(self, S, K, r, q, T, sigma)
        self.price = self.K * math.exp(-self.r * self.T) * norm.cdf(-self.d2) - self.S * math.exp(-self.q * self.T) * norm.cdf(-self.d1)
        self.get_greeks()
            
    def get_greeks(self):
        self.delta = - math.exp(-self.q * self.T) * norm.cdf(-self.d1)
        self.gamma = (math.exp(-self.q * self.T) * norm.pdf(self.d1)) / (self.S * self.sigma * math.sqrt(self.T))
        self.vega = self.S * math.exp(-self.q * self.T) * norm.pdf(self.d1) * math.sqrt(self.T)
        self.theta = (- math.exp(-self.q * self.T) * self.S * norm.pdf(self.d1) * self.sigma / (2 * math.sqrt(self.T))) + self.r * self.K * math.exp(-self.r * self.T) * norm.cdf(-self.d2) - self.q * self.S * math.exp(-self.q * self.T) * norm.cdf(-self.d1)
        
class BSEuropeanCall(BSOption):
    """ BS call options"""
    def __init__(self, S, K, r, q, T, sigma):
        BSOption.__init__(self, S, K, r, q, T, sigma)
        self.price = self.S * math.exp(-self.q * self.T) * norm.cdf(self.d1) - self.K * math.exp(-self.r * self.T) * norm.cdf(self.d2)
        self.get_greeks()
            
    def get_greeks(self):
        self.delta = math.exp(-self.q * self.T) * norm.cdf(self.d1)
        self.gamma = (math.exp(-self.q * self.T) * norm.pdf(self.d1)) / (self.S * self.sigma * math.sqrt(self.T))
        self.vega = self.S * math.exp(-self.q * self.T) * norm.pdf(self.d1) * math.sqrt(self.T)
        self.theta = (- math.exp(-self.q * self.T) * self.S * norm.pdf(self.d1) * self.sigma / (2 * math.sqrt(self.T))) - self.r * self.K * math.exp(-self.r * self.T) * norm.cdf(self.d2) + self.q * self.S * math.exp(-self.q * self.T) * norm.cdf(self.d1)
    
        
class DataLoader:
    """ Class used to load mkt data from excel"""
    
    def __init__(self):
        self.vix = False
        self.sofr = False
        self.IVV = False
        
        self.get_vix_data()
        self.get_sofr_data()
        self.get_underlying()
        
        self.dates = [date for date in self.IVV["asofdate"] if ((date in list(self.sofr["asofdate"])) & (date in list(self.vix.asofdate)))]
        self.dates.sort()
        
        self.sofr = self.sofr.loc[self.sofr["asofdate"].isin(self.dates)].copy()
        self.vix = self.vix.loc[self.vix["asofdate"].isin(self.dates)].copy()
        self.IVV = self.IVV.loc[self.IVV["asofdate"].isin(self.dates)].copy()

        
    def get_vix_data(self):
        self.vix = pd.read_excel(os.path.join(os.getcwd(), "data", 'vix.xlsx') , index_col = 0)
        self.vix["asofdate"] = (pd.to_datetime(self.vix["asofdate"])).dt.date
        self.vix = self.vix.loc[self.vix.asofdate >= datetime.date(2020, 1, 2)].copy()
        self.vix.reset_index(inplace = True, drop = True)

    def get_sofr_data(self):
        self.sofr = pd.read_excel(os.path.join(os.getcwd(), "data", 'sofr.xlsx'))
        self.sofr["asofdate"] = (pd.to_datetime(self.sofr["Effective Date"])).dt.date
        self.sofr = self.sofr.loc[self.sofr["Rate Type"] == "SOFR"].copy()
        self.sofr.reset_index(inplace = True, drop = True)
        
    def get_underlying(self):
        self.IVV = yf.download(["IVV"] , start = datetime.date(2020, 1, 2), end=datetime.date(2023, 6, 1), group_by='tickers')
        self.IVV["asofdate"] = (pd.to_datetime(self.IVV.index)).date

    def get_vix_contract(self, vix_, trade_date_, time_to_expiration_):
        row_val = vix_.loc[(vix_.asofdate == trade_date_)].apply(lambda x: abs(x.YTTE- time_to_expiration_), axis=1).idxmin()
        return float(vix_.loc[row_val, "Settle"])/100
