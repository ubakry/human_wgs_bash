"""
Update Workflow Status Script 
Copyright 2021 Usama Bakry (u.bakry@icloud.com)

This file is used for updating workflow status in database.
"""

# Importing libraries
import argparse
from os import sep
import sys
from datetime import datetime
import pandas as pd
import pymysql
from sqlalchemy import create_engine

args = None
time_fmt = '%H:%M:%S'

# ----------------------------------------------------------------------
def get_args():
    
    """
	get_args() function creates required and optional arguments and returns a copy of the argument list to be used within the script.
	"""
    parser = argparse.ArgumentParser(
        description="This script is used for updating workflow status in database.",
        epilog="This is where you might put example usage"
        )

    # sample args
    parser.add_argument('-p', action="store", required=True, help='Input project name')
    parser.add_argument('-r', action="store", required=True, help='Input run name')
    parser.add_argument('-s', action="store", required=True, help='Input sample name')

    # workflow args
    parser.add_argument('--convert', action="store", default='In Progress' ,help="Convert BCL to fastq status (Waiting, In Progress, Fail, or [duration] - default=\'In Progress\')")
    parser.add_argument('--qc', action="store", default='Waiting' ,help='Quality control status (Waiting, In Progress, Fail, or [duration]) - default=\'Waiting\'')
    parser.add_argument('--alignment', action="store", default='Waiting' ,help='Alignment status (Waiting, In Progress, Fail, or [duration]) - default=\'Waiting\'')
    parser.add_argument('--mark_duplicates', action="store", default='Waiting' ,help='Mark duplicates status (Waiting, In Progress, Fail, or [duration]) - default=\'Waiting\'')
    parser.add_argument('--variant_calling', action="store", default='Waiting' ,help='Variant calling status (Waiting, In Progress, Fail, or [duration]) - default=\'Waiting\'')
    parser.add_argument('--annotation', action="store", default='Waiting' ,help='Annotation status (Waiting, In Progress, Fail, or [duration]) - default=\'Waiting\'')
    parser.add_argument('--prs_calculations', action="store", default='Waiting' ,help='PRS calculations status (Waiting, In Progress, Fail, or [duration]) - default=\'Waiting\'')
    parser.add_argument('--report_generation', action="store", default='Waiting' ,help='Report generation status (Waiting, In Progress, Fail, or [duration]) - default=\'Waiting\'')
    
    arguments = vars(parser.parse_args())
    return arguments
# ----------------------------------------------------------------------

def get_current_timestamp():
    """
	Get current time and format it
	"""
    ts_now = datetime.now()
    ts_format = str(ts_now.day)+"/"+str(ts_now.month)+"/"+str(ts_now.year)+" "+str(ts_now.hour)+":"+str(ts_now.minute)+":"+str(ts_now.second)
    return ts_format

# ----------------------------------------------------------------------
def main():
    """
	The main script
	"""
    # Database connection info
    table_name= 'workflow_status'
    sql_engine= create_engine('mysql+pymysql://root:QSZWAX191994ub$@localhost:3306/bioinfo_operations')
    connection    = sql_engine.connect()

    # Get the table into data frame
    try:
        get_df = pd.read_sql(table_name,connection,columns=['project_id','run_id','sample_name'])
    except ValueError as vx:
        print("Value Error:", vx)
    except Exception as ex:
        print("Exception:", ex)
    else:
        print("Getting table is done successfully.")

    # Check if entry with input project_id, run_id and sample_name is exist
    if (get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s'])].shape[0] == 0):

        # Convert arguments to data frame
        row_data = {
            'project_id' : [args['p']],
            'run_id' : [args['r']],
            'sample_name' : [args['s']],
            'start_timestamp' : [get_current_timestamp()],
            'convert_bcl' : [args['convert']],
            'qc' : [args['qc']],
            'alignment' : [args['alignment']],
            'mark_duplicates' : [args['mark_duplicates']],
            'variant_calling' : [args['variant_calling']],
            'annotation' : [args['annotation']],
            'prs_calculations' : [args['prs_calculations']],
            'report_generation' : [args['report_generation']],
            'end_timestamp' : ["In Progress"],
            'duration' : ["In Progress"]
        }

        row_df = pd.DataFrame(data=row_data)

        # Add row to database
        try:
            frame = row_df.to_sql(table_name,connection, if_exists='append', index= False)
        except ValueError as vx:
            print("Value Error:", vx)
        except Exception as ex:
            print("Exception:", ex)
        else:
            print("Row is added successfully.")
        finally:
            connection.close()

    else:
        # Get the table into data frame
        try:
            get_df = pd.read_sql(table_name,connection)
        except ValueError as vx:
            print("Value Error:", vx)
        except Exception as ex:
            print("Exception:", ex)
        else:
            print("Getting table is done successfully.")

        # selected_row = get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s'])]

        # Check values to update data
        if (args['convert'] != "In Progress"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['convert_bcl']] = args['convert']

        if (args['qc'] != "Waiting"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['qc']] = args['qc']

        if (args['alignment'] != "Waiting"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['alignment']] = args['alignment']

        if (args['mark_duplicates'] != "Waiting"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['mark_duplicates']] = args['mark_duplicates']

        if (args['variant_calling'] != "Waiting"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['variant_calling']] = args['variant_calling']

        if (args['annotation'] != "Waiting"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['annotation']] = args['annotation']

        if (args['prs_calculations'] != "Waiting"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['prs_calculations']] = args['prs_calculations']

        if (args['report_generation'] != "Waiting"):
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['report_generation']] = args['report_generation']

        selected_row = get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s'])]

        # End timestamp
        if (selected_row.iloc[0,11] == "Waiting" or selected_row.iloc[0,11] == "In Progress"):
            end_timestamp = 'In Progress'  
        else:
            end_timestamp = get_current_timestamp()
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['end_timestamp']] = end_timestamp

        # Duration
        if (end_timestamp == 'In Progress'):
            duration = 'In Progress'
        else:
            start_datetime_ls = selected_row.iloc[0,3].split(" ")[0].split("/") + selected_row.iloc[0,3].split(" ")[1].split(":")
            start_datetime_dt = datetime(int(start_datetime_ls[2]),int(start_datetime_ls[1]),int(start_datetime_ls[0]),int(start_datetime_ls[3]),int(start_datetime_ls[4]),int(start_datetime_ls[5]))
            
            end_datetime_ls = end_timestamp.split(" ")[0].split("/") + end_timestamp.split(" ")[1].split(":")
            end_datetime_dt = datetime(int(end_datetime_ls[2]),int(end_datetime_ls[1]),int(end_datetime_ls[0]),int(end_datetime_ls[3]),int(end_datetime_ls[4]),int(end_datetime_ls[5]))

            tdelta = end_datetime_dt - start_datetime_dt
            duration = str(int(round(tdelta.total_seconds()/60,0))) +"m"
            get_df.loc[(get_df['project_id'] == args['p']) & (get_df['run_id'] == args['r']) & (get_df['sample_name'] == args['s']),['duration']] = duration
            

        # Add row to database
        try:
            frame = get_df.to_sql(table_name,connection, if_exists='replace', index= False)
        except ValueError as vx:
            print("Value Error:", vx)
        except Exception as ex:
            print("Exception:", ex)
        else:
            print("Update is successful.")
        finally:
            connection.close()
# ----------------------------------------------------------------------

if __name__ == '__main__':
	args = get_args()
	main()

