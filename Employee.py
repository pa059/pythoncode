class Employee:    
     # init method initializes the name and id attributes
    def __init__(self, name, idNumber):
        self.__name = name
        self.__idNumber = idNumber

    # The set method sets the name attribute
    def set_name(self, name):
        self.__name = name

    # The set method sets the id number attribute
    def set_idNumber(self, idNumber):
        self.__idNumber = idNumber

    # The get method returns the name attribute
    def get_name(self):
        return self.__name

    # The get method returns the id number attribute
    def get_idNumber(self):
        return self.__idNumber

# The ProductionWorker class is a subclass of the Employee class
class ProductionWorker(Employee):

    # init method initializes all the attributes and inherits the super class' attributes
    def __init__(self, name, idNumber, shiftNumber, hourlyPayRate):
        Employee.__init__(self, name, idNumber)
        self.__shiftNumber = shiftNumber
        self.__hourlyPayRate = hourlyPayRate

    # The set method sets the shift number attribute
    def set_shiftNumber(self, shiftNumber):
        self.__shiftNumber = shiftNumber

    # The set method sets the hourly pay attribute
    def set_hourlyPayRate(self, hourlyPayRate):
        self.__hourlyPayRate = hourlyPayRate,'.2f'

    # The get method returns the shift attribute
    def get_shiftNumber(self):

        # if user selects 1 the attibute is converted into a string when displayed
        if self.__shiftNumber == 1:
            self.__shiftNumber = 'Day shift'
        elif self.__shiftNumber == 2:
            self.__shiftNumber = 'Night shift'
        else:
            print('That is not a valid shift')
        return self.__shiftNumber

    # The get method returns the hourly pay rate attribute
    def get_hourlyPayRate(self):
        return self.__hourlyPayRate

    # The string method converts numerical data into a string and returns it
    def __str__(self):
        result = "Name: " + self.get_name() + "\n" + \
                "ID number: " + str(self.get_idNumber()) + \
                "\nShift: " + self.get_shiftNumber() + \
                "Hourly pay: " + str(self.get_hourlyPayRate())
        return result

# The shift supervisor class is a subclass of the Employee class
class ShiftSupervisor(Employee):

    # init method initializes attributes and inherits the employee class' attributes
    def __init__(self, name, idNumber, salary, bonus):
        Employee.__init__(self, name, idNumber)
        self.__salary = salary
        self.__bonus = bonus

    # The set method sets the salary attribute
    def set_salary(self, salary):
        self.__salary = salary

    # The set method sets the bonus attribute
    def set_bonus(self, bonus):
        self.__bonus = bonus

    # The get method returns the salary attribute
    def get_salary(self):
        return self.__salary

    # The get method returns the bonus attribute
    def get_bonus(self):
        return self.__bonus