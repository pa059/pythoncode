# imports employee class
import Employee

# Main function
def main():

    # Calls the function to make worker list
    workers = makeWorkerList()

    # Output
    print('Here are all the workers hired today:')
    display_worker(workers)

# Function to take user input, store in the object, then store in a list
def makeWorkerList():

    # Creates an empty list to store employee data
    worker_list = []

    # Input to get the number of worker's data the user will be entering in
    numOfWorkers = int(input('How many workers hired today? '))

    # Loop to iterate through each worker
    for num in range(numOfWorkers):
        print('Enter the data for worker', int(num + 1))

        # Input
        name = input('Enter the name of the worker: ')
        idNumber = int(input('Enter the id number: '))
        shiftNumber = int(input('Enter the shift number (1 or 2): '))
        hourlyPayRate = float(input('Enter the hourly pay: '))
        print()

        # Creates the Production worker object to store the atributes in
        worker = Employee.ProductionWorker(name, idNumber, shiftNumber, hourlyPayRate)

        # append the data into a list
        worker_list.append(worker)

    return worker_list

# function to display the information
def display_worker(worker_list):
    for worker in worker_list:
        print('Name:', worker.get_name())
        print('ID number:', worker.get_idNumber())
        print('Shift:',worker.get_shiftNumber())
        print('Hourly rate: $',worker.get_hourlyPayRate())
        print()

# Calls the main function
main()