// Class to represent the environment
class Verden { // World
    boolean regner; // raining
    int dag; // 1 is Monday, 2 is Tuesday, ..., 7 is Sunday, 8 is Monday...

    // Method to check if it's Sunday
    public boolean erSondag() { // isSunday
        return dag % 7 == 0; // Returns true if the day is Sunday
    }
}

// Class to represent the Robot
class Robot {
    private String navn; // name
    private double batteriNivaa; // batteryLevel
    private int avstandTilParken; // distanceToPark in meters
    private String botType; // botType

    // Constructor to initialize Robot properties
    public Robot(String navn, double batteriNivaa, int avstandTilParken, String botType) {
        this.navn = navn; // name
        setBatteriNivaa(batteriNivaa); // Ensures valid battery level
        this.avstandTilParken = avstandTilParken; // distanceToPark
        this.botType = botType; // botType
    }

    // Getters and Setters
    public String getNavn() { // getName
        return navn;
    }

    public void setNavn(String navn) { // setName
        this.navn = navn;
    }

    public double getBatteriNivaa() { // getBatteryLevel
        return batteriNivaa;
    }

    public void setBatteriNivaa(double batteriNivaa) { // setBatteryLevel
        // Validates battery level between 0 and 100
        this.batteriNivaa = Math.max(0, Math.min(batteriNivaa, 100));
    }

    public int getAvstandTilParken() { // getDistanceToPark
        return avstandTilParken;
    }

    public void setAvstandTilParken(int avstandTilParken) { // setDistanceToPark
        this.avstandTilParken = Math.max(0, avstandTilParken); // Ensures distance is not negative
    }

    public String getBotType() { // getBotType
        return botType;
    }

    public void setBotType(String botType) { // setBotType
        this.botType = botType;
    }

    // Method to return the robot's status
    public String giStatus() { // giveStatus
        return "Dette er bot " + navn + " av type " + botType + ". Denne enheten har " + batteriNivaa +
               "% batterikapasitet igjen og bor " + avstandTilParken + " meter fra parken."; 
               // This is bot {name} of type {botType}. This unit has {batteryLevel}% battery capacity left 
               // and lives {distanceToPark} meters from the park.
    }

    // Method to check if the robot can go to the park based on the environment
    public boolean gaaTilParken(Verden verden) { // goToPark(World world)
        System.out.println(navn + " sjekker om det er mulig å gå til parken...");
        // {name} is checking if it is possible to go to the park...

        // Check if it's raining
        if (verden.regner) { // world.isRaining
            System.out.println("Det regner, så " + navn + " kan ikke gå til parken.");
            // It's raining, so {name} cannot go to the park.
            return false;
        }

        // Check if it's Sunday. Can only go to the park on Sundays.
        if (!verden.erSondag()) { // world.isSunday
            System.out.println("Det er ikke søndag, så " + navn + " kan ikke gå til parken.");
            // It is not Sunday, so {name} cannot go to the park.
            return false;
        }

        // Check battery level. A bot needs 1% battery per 100 meters to go to the park.
        double neededBattery = avstandTilParken / 100.0; // distanceToPark / 100.0
        if (batteriNivaa < neededBattery) { // batteryLevel < neededBattery
            System.out.println(navn + " har ikke nok batteri til å gå til parken.");
            // {name} does not have enough battery to go to the park.
            return false;
        }

        // Reduce battery after going to the park
        setBatteriNivaa(batteriNivaa - neededBattery); // setBatteryLevel
        System.out.println(navn + " har gått til parken og nå har " + batteriNivaa + "% batterikapasitet igjen.");
        // {name} has gone to the park and now has {batteryLevel}% battery capacity left.
        return true;
    }

    // Method to check if the robot can attend the dance club
    public boolean deltaIDanseklubb(Verden verden) { // joinDanceClub(World world)
        // Bot type must be "B-Bot" and need a minimum of 50% battery
        // Dance club is closed on Mondays (day 1)
        if (botType.equals("B-Bot") && batteriNivaa >= 50 && verden.dag % 7 != 1) {
            System.out.println(navn + " kan delta i danseklubben.");
            // {name} can join the dance club.
            return true;
        } else {
            System.out.println(navn + " kan ikke delta i danseklubben.");
            // {name} cannot join the dance club.
            return false;
        }
    }
}

// Main class to test the Robot and Verden classes
public class RobotOppgave { // RobotTask
    public static void main(String[] args) {
        // Creating an environment object (Verden)
        Verden dagensVerden = new Verden(); // today'sWorld
        dagensVerden.regner = false; // it's not raining
        dagensVerden.dag = 7; // For example, Sunday

        // Creating two Robot objects
        Robot robot1 = new Robot("Dancatron 4000", 80.0, 1500, "B-Bot");
        Robot robot2 = new Robot("Spark-E", 30.0, 500, "Toaster");

        // Print the status of both robots
        System.out.println(robot1.giStatus()); // giveStatus
        System.out.println(robot2.giStatus());

        // Check if the robots can attend the dance club
        robot1.deltaIDanseklubb(dagensVerden); // joinDanceClub
        robot2.deltaIDanseklubb(dagensVerden);

        // Check if the robots can go to the park
        robot1.gaaTilParken(dagensVerden); // goToPark
        robot2.gaaTilParken(dagensVerden);
    }
}
