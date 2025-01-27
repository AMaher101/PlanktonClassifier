import os
from datetime import datetime
from classifier import Classifier

csvs = sorted(os.listdir(f"{os.getcwd()}/inputs"))

for csv_name in csvs:
    # classify mixotrophs
    classified = Classifier(csv_name)
    
    # save dataframe to excel
    classified.all_classified.to_excel(f"outputs/{csv_name}-{str(datetime.now())}.xlsx")
    print(csv_name + " done.")