
import Employee

def main():

    supervisor= getSupervisor()
    displaySupervisor(supervisor)


def getSupervisor():

    name=input("Enter the supervisor name:")
    idNumber=int(input("Enter the supervisor id:"))
    salary = float(input("Enter the salary:"))
    bonus=float(input("Enter the bonus:"))

    supervisor = Employee.ShiftSupervisor(name,idNumber,salary,bonus)
    return supervisor

def displaySupervisor(supervisor):

    print("Shift supervisor")
    print("************")
    print("Name:", supervisor.get_name())
    print("ID:", supervisor.get_idNumber())
    print("Salary",supervisor.get_salary())
    print("Bonus",supervisor.get_bonus())

main()
