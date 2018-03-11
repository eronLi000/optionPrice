import pandas_datareader.data as pdr

source='yahoo'
def getStockPrice(valueDate,symbol=""):
    if symbol =="":
        self.msg = QMessageBox.question(self,'Invalid Stock',"Sorry, please select a stock.", QMessageBox.Ok)
        return
    stock = pdr.DataReader(symbol, source,
                            valueDate,
                            valueDate)
    stockOpenPrice = round(stock.ix[str(valueDate)]["Open"],2)
    return stockOpenPrice