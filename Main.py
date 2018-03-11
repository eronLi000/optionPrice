import sys
from PyQt5.QtWidgets import *
import PyQt5
#import math
import StockPrice
import FileReader
import datetime
import time
from PyQt5 import QtCore, uic
import numpy as np
from scipy.stats import norm 
from datetime import date

QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
ui_MW, QtBaseClass = uic.loadUiType("Main.ui")
ui_SL, QtBaseClass = uic.loadUiType("StockList.ui")

#Main Window
class Ui_MainWindow(QMainWindow, ui_MW):
    def __init__(self, parent=None):
        super(Ui_MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.search_button.clicked.connect(self.popStockList)
        self.calculate_button.clicked.connect(self.calculatePrice)
        self.reset_button.clicked.connect(self.resetValue)
		
    def getValueDate(self):
        return datetime.date(self.value_date.date().year(),self.value_date.date().month(),self.value_date.date().day())
    
    def resetValue(self):
        self.stock_price.setText("")
        self.volatility.setText("")
        self.strike_price.setText("")
        self.interest_rate.setText("")
        self.dividend_yield.setText("")
        self.stock_price.setText("")
        self.M_value.setText("")
        self.N_value.setText("")
        self.theoretical_value.setText("")

    def calculatePrice(self):
        #get stock price
        S,K,V,R,Q,T,M,N=0.0,0.0,0.0,0.0,0.0,0.0,0,0
        try:
            S=float(self.stock_price.text())
        except:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for stock price is invalid.", QMessageBox.Ok)
        #get stock price        
        try:
            K=float(self.strike_price.text())
        except:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for strike price is invalid.", QMessageBox.Ok)
        #get volatility rate        
        try:
            V=float(self.volatility.text())/100
        except:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for volatility is invalid.", QMessageBox.Ok)
        #get interest rate        
        try:
           R=(float(self.interest_rate.text()))/100
        except:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for interest rate is invalid.", QMessageBox.Ok)
        #get dividend yield        
        try:
            Q=(float(self.dividend_yield.text()))/100
        except:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for dividend_yield is invalid.", QMessageBox.Ok)
        #get M value        
        try:
            M=int(self.M_value.text())
        except:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for M is invalid.", QMessageBox.Ok)
        #get N value        
        try:
            N=int(self.N_value.text())
        except:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for N is invalid.", QMessageBox.Ok)

        #get operationType
        operationType = self.operation_type.currentText()
        #get pde solution type
        pdeType = self.pde_solution.currentText()
        #get T
        dateValue = datetime.date(self.value_date.date().year(),self.value_date.date().month(),self.value_date.date().day())
        expireDate = datetime.date(self.expiration_date.date().year(),self.expiration_date.date().month(),self.expiration_date.date().day()) 
        if expireDate < dateValue:
            self.msg = QMessageBox.question(self,'Invalid Input',"Sorry, Input for expiration_date is invalid.", QMessageBox.Ok)
        T=(expireDate-dateValue).days/365
        

        result = 0
        d1=(np.log(S/K)+(R-Q+V**2/2)*T)/(V*np.sqrt(T))
        d2 = d1-(V*np.sqrt(T))
        Smax = 2*K
        dS = Smax/M
        dt = T/N
        C=round(S*norm.cdf(d1)*np.exp(-Q*T)-norm.cdf(d2)*K*np.exp(-R*T),2)
        p=round(K*np.exp(-R*T)*norm.cdf(-d2)-S*np.exp(-Q*T)*norm.cdf(-d1),2)
        k=int(S/dS)
        b= np.zeros((M+1,M+1))
        
        for x in range(1,M):
            b[x][x] = (1)-dt*(V**2*x**2+R)
            b[x][x-1] = (1/2)*dt*(V**2*x**2-R*x)
            b[x][x+1] = (1/2)*dt*(V**2*x**2+R*x)
        b[0][0] = 1
        b[M][M] = 1
        np_b = np.array(b)
        
        c = np.zeros((M+1,M+1))
        for x in range(1,M):
            c[x][x]= (1)+dt*(V**2*x**2+R)
            c[x][x-1]= (1/2)*dt*((R-Q)*x-V**2*x**2)
            c[x][x+1]= (-1/2)*dt*((R-Q)*x+V**2*x**2)
        c[0][0]=1
        c[M][M]=1
        np_c=np.array(c)
        np_cinv= np.linalg.inv(np_c)

        d = np.zeros((M+1,M+1))
        for x in range(1,M):
            d[x][x]= 1+(-1/2)*dt*(V**2*x**2+R)
            d[x][x-1]= (1/4)*dt*(V**2*x**2-(R-Q)*x)
            d[x][x+1]= (1/4)*dt*((R-Q)*x+V**2*x**2)
        d[0][0]=1
        d[M][M]=1
        np_d=np.array(d)

        e = np.zeros((M+1,M+1))
        for x in range(1,M):
            e[x][x]= 1-(-1/2)*dt*(V**2*x**2+R)
            e[x][x-1]= (-1/4)*dt*(V**2*x**2-(R-Q)*x)
            e[x][x+1]= (-1/4)*dt*((R-Q)*x+V**2*x**2)
        e[0][0]=1
        e[M][M]=1
        np_e=np.array(e)
        np_einv= np.linalg.inv(np_e)

        table = np.zeros((M+1,N+1))
        tableC = np.array(table)
        tableP = np.array(table)
        fmax=np.array([i*dS for i in range(M+1)])
        tableC[:,N] = np.maximum(fmax-K,0)
        tableP[:,N] = np.maximum(K-fmax,0)
        np_tableEC=np.array(tableC)
        np_tableEP=np.array(tableP)
        np_tableIC=np.array(tableC)
        np_tableIP=np.array(tableP)
        np_tableCNC=np.array(tableC)
        np_tableCNP=np.array(tableP)

        for y in range(N-1,-1,-1):
            np_tableEC[:,y]=np_b.dot(np_tableEC[:,y+1])
            np_tableEP[:,y]=np_b.dot(np_tableEP[:,y+1])
            np_tableIC[:,y]=(np_cinv).dot(np_tableIC[:,y+1])
            np_tableIP[:,y]=(np_cinv).dot(np_tableIP[:,y+1])
            vbc=(np_d).dot(np_tableCNC[:,y+1])
            np_tableCNC[:,y]=np_einv.dot(vbc)
            vbp=(np_d).dot(np_tableCNP[:,y+1])
            np_tableCNP[:,y]=np_einv.dot(vbp)

        if operationType == "Call Operation" and pdeType=="Crank-nicolson":
            COCN = round( np_tableCNC[k,0]+(np_tableCNC[k+1,0]-
                np_tableCNC[k,0])/dS*(S-k*dS),2)
            self.theoretical_value.setText(str(COCN))
        elif operationType == "Call Operation" and pdeType=="Implicit":
            ICO =round(np_tableIC[k,0]+(np_tableIC[k+1,0]-
                np_tableIC[k,0])/dS*(S-k*dS),2)        
            self.theoretical_value.setText(str(ICO))
        elif operationType == "Call Operation" and pdeType=="Explicit":
            ECO =round(np_tableEC[k,0]+(np_tableEC[k+1,0]-
                np_tableEC[k,0])/dS*(S-k*dS),2)
            self.theoretical_value.setText(str(ECO))
        elif operationType == "Pull Operation" and pdeType=="Crank-nicolson":
            CNPO = round( np_tableCNP[k,0]+(np_tableCNP[k+1,0]-
                np_tableCNP[k,0])/dS*(S-k*dS),2)
            self.theoretical_value.setText(str(CNPO))
        elif operationType == "Pull Operation" and pdeType=="Implicit":
            IPO =round(np_tableIP[k,0]+(np_tableIP[k+1,0]-
                np_tableIP[k,0])/dS*(S-k*dS),2)
            self.theoretical_value.setText(str(IPO))
        else:#if operationType == "Pull Operation" and pdeType=="Explicit":
            EPO =round(np_tableEP[k,0]+(np_tableEP[k+1,0]-
                np_tableEP[k,0])/dS*(S-k*dS),2)
            self.theoretical_value.setText(str(EPO))

    def popStockList(self):
        self.dialog = None
        try:
            self.dialog = Ui_StockDialog(self,ui_SL)
            self.dialog.show()
        except ApplicationException:
            self.msg = QMessageBox.question(self,'No Stock Found',"Sorry, we couldn't find what you were looking for.", QMessageBox.Ok)

class Ui_StockDialog(QDialog,ui_SL):
    def __init__(self, parent=None, data=None):
        super(Ui_StockDialog, self).__init__(parent)
        self.setupUi(self)
        self.setTableData()
        self.show()
        
        self.okButton.clicked.connect(self.okBtn)
        self.cancelButton.clicked.connect(self.closeDialog)
       
    def setTableData(self):
        stocklist = FileReader.scan_file('stocklist.csv')
        if len(stocklist) == 0:
            self.msg = QMessageBox.question(self, 'Not enough stock', "no stock.",QMessageBox.Ok)
        index = 0
        for key, value in stocklist.items():
            self.stocklistTable.insertRow(index)
            self.stocklistTable.setItem(index, 0, QTableWidgetItem(key.strip("\"")))
            self.stocklistTable.setItem(index, 1, QTableWidgetItem(value.strip("\"")))
            index+=1

    def closeDialog(self):
        return self.accept()
    
    def okBtn(self):
        if self.stocklistTable.selectedIndexes() == []:
            self.msg = QMessageBox.question(self, 'Alert',"Please select stock", QMessageBox.Ok)
            return
        index = self.stocklistTable.selectedIndexes()
        #print("index:",index)
        symbol = index[0].data()
        stockname = index[1].data()
        main = self.parent()
        main.stock_name.setText(symbol)
        dateValue = main.getValueDate()
        try:
            stockPrice = StockPrice.getStockPrice(dateValue,symbol)
            main.stock_price.setText(str(stockPrice))
        except:
            self.msg = QMessageBox.question(self, 'Price unavailable',"Price unavailable on the selected day", QMessageBox.Ok)
        self.closeDialog()

def main():
    app = QApplication(sys.argv)
    main = Ui_MainWindow()
    main.show()
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    main()

    '''
        if pdeType == "Crank-nicolson":
            result = bs_crank_nic(S,K,r,sigma,T,M,N,operationType)
        elif pdeType == "Implicit":
            result = bs_implicit(S,K,r,sigma,T,M,N,operationType)
        else:
            result = bs_explicit(S,K,r,sigma,T,M,N,operationType)
        self.theoretical_value.setText(str(result))

    def bs_explicit(S,K,r,sigma,T,M,N,operationType):
        dlata_T = T/N
        Smax = 2*K
        delta_S = Smax/M

        if operationType == "Call Operation":
            for i 
        else:
        
    def bs_crank_nic(S,K,r,sigma,T,M,N,operationType):
    def bs_implicit(S,K,r,sigma,T,M,N,operationType):
    '''