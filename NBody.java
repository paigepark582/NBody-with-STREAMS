import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;
import java.util.stream.*;
import java.io.*;
import java.util.concurrent.*;
import java.lang.Math;
import javax.swing.*;
import javax.swing.Timer;
import java.awt.geom.*;

public class NBody extends JFrame {
    public static double G = 6.73e-4;
    public static double DT = .5;
    public static double mass = 100.0;
    public static final int framesize = 1000;

    public static int numWorkers;
    public static Semaphore[][] stages;
    public static Semaphore forMain;
    public static Semaphore forWorker;
    public static int numStages;
    public static long collisions;
    public static Object collisionSync = 0;

    public static volatile List<Body> bodies;
    public NBody() { //nbody.repaint()
        setSize(framesize, framesize);
        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        setVisible(true);
        add(new JPanel() {
            @Override
            public void paintComponent(final Graphics graphics){
                super.paintComponent(graphics);
                graphics.setFont(new Font(graphics.getFont().getFontName(), Font.PLAIN, 11));
                Graphics2D g2 = (Graphics2D) graphics;
                NBody.bodies.stream().forEach(b -> {
                    g2.draw(new Ellipse2D.Double(b.getX() - (b.radius),
                            b.getY() - (b.radius), 2.0 * b.radius, 2.0 * b.radius));
                    String[] str=b.label.split("\\n");
                    float x = (float) b.getX();
                    float y = (float) b.getY();
                    g2.drawString(str[0], (float) (x - b.radius), y - (float)b.radius);
                    g2.drawString(str[1], (float) (x - b.radius), (float) (y + b.radius * 4.0/3.0 - b.radius));
                    g2.drawString(str[2], (float) (x - b.radius), (float) (y + b.radius * 2.0/3.0 - b.radius));
                });
            }
        });
        new Timer(1000 / 60, new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                repaint();
            }
        }).start();
    }
    public static Random gen;
    public static void main(String[] args){
        if(args.length < 4){
            System.out.println("Usage NBody <workers> <bodies> <radius> <time steps>");
            System.out.println("Use --help to see all options (must come after 4 args)");
            System.exit(1);
        }
        if(Arrays.stream(args).anyMatch("--help"::equals)) {
            System.out.println("-seed <num> : a seed for the random generator");
            System.out.println("-radRange <low>-<high> : randomize the radius of each circle (adds a number in the range to the 3rd argument)");
            System.out.println(String.format("-massRange <low>-<high> : randomize the mass of each circle (default: %.2f)",mass));
            System.out.println(String.format("-DT <interval> : time interval (default: %.2f)",DT));
            System.out.println(String.format("-g : turn on graphics"));
            System.exit(1);
        }
        numWorkers = Integer.parseInt(args[0]);
        int numBodies = Integer.parseInt(args[1]);
        double size = Double.parseDouble(args[2]);
        int steps = Integer.parseInt(args[3]);
        bodies = new ArrayList<Body>();
        double x=0, y=0;
        int rlow=0, rhigh=1, mlow=0, mhigh=1;
        gen = new Random();
        JFrame frame = null;
        if(args.length > 4){
            for(int i = 4; i < args.length; i++) {
                String[] temp;
                switch (args[i]) {
                    case "-seed":
                        i++;
                        gen = new Random(Integer.parseInt(args[i]));
                        break;
                    case "-radRange":
                        i++;
                        temp = args[i].split("-");
                        rlow = Integer.parseInt(temp[0]);
                        rhigh = Integer.parseInt(temp[1]);
                        break;
                    case "-massRange":
                        i++;
                        temp = args[i].split("-");
                        mlow = Integer.parseInt(temp[0]);
                        mhigh = Integer.parseInt(temp[1]);
                        break;
                    case "-DT":
                        i++;
                        DT = Double.parseDouble(args[i]);
                        break;
                    case "-g":
                        frame = new NBody();
                        break;
                    default:
                        System.out.println("Unrecognized Option: " + args[i]);
                }
            }
        }
        for(int i = 0; i < numBodies; i++){
            double rad = size + gen.nextInt(rhigh - rlow);
            double m = mass + gen.nextInt(mhigh - mlow);
            Body b = new Body(gen.nextDouble() * (framesize - 2 * rad) + rad, gen.nextDouble() * (framesize - 2 * rad) + rad, rad, m);
            //no wall collisions to start/collisions in general
//            System.out.println(String.format("%.2f, %.2f",b.getX(), b.getY()));
            if (bodies.stream().anyMatch( body -> hasCollided(b, body))) {
                i -= 1;
                continue;
            }
            bodies.add(b);
        }

        numStages = (int)Math.ceil(Math.log(numWorkers)/Math.log(2));
        stages = new Semaphore[numWorkers][numStages];
        for(int i = 0; i < numWorkers; i++) {
            for(int j = 0; j < numStages; j++) {
                stages[i][j] = new Semaphore(0);
            }
        }
        forMain = new Semaphore(0);
        forWorker = new Semaphore(0);
        
        long startTime = System.nanoTime(), endTime; 

        NBodyWorker[] workers = new NBodyWorker[numWorkers];
        int split = numBodies/numWorkers;
        for(int i = 0; i < numWorkers; i++) {
            workers[i] = new NBodyWorker(i, i*split, (i != numWorkers - 1) ? ((i+1)*split) : (numBodies), steps);
            workers[i].start();
        }
        endTime=System.nanoTime();
        long micros=(endTime-startTime)/1000;
        int seconds=0;
        while(micros>1000000){
            micros-=1000000;
            seconds++;
        }
        System.out.println(String.format("computation time: %d seconds, %d microseconds (for creating processes)",seconds,micros));
        
        startTime = System.nanoTime();
        try {
            for (int i = 0; i < steps; i++) {
                forMain.acquire();
                bodies = Collections.unmodifiableList(
                        Arrays.stream(workers)
                        .flatMap(w -> w.bodies.stream())
                        .collect(Collectors.toList()));
                //System.out.println(bodies.size());
                forWorker.release();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        endTime=System.nanoTime();
        micros=(endTime-startTime)/1000;
        seconds=0;
        while(micros>1000000){
            micros-=1000000;
            seconds++;
        }
        System.out.println(String.format("computation time: %d seconds, %d microseconds",seconds,micros));
        System.out.println(String.format("Collisions detected: %d",collisions/2));
        try {
            File f = new File("Results.txt");
            FileWriter out = new FileWriter(f, false);
            out.write("Positions and Velocities\n");
            bodies.stream()
                .forEach(b -> {
                    try{out.write(String.format("(%.2f, %.2f) vx: %.2f vy: %.2f\n", b.getX(), b.getY(), b.vx, b.vy));} 
                    catch (IOException e) {
                        e.printStackTrace();
                        System.exit(1);}});
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        if(frame != null) {
            frame.setVisible(false);
            frame.dispose();
        }
        System.exit(1);
    }
    public static ArrayList<Body> calculateForces(List<Body> bodies) {
        return bodies.stream()
                .map(NBody::calculateAgainst)
                .collect(Collectors.toCollection(ArrayList::new));
    }
    private static Body calculateAgainst(Body b) { //bodies
        Collector<Body, Accum, Body> forceCollector =
                Collector.of(
                        () -> new Accum(0.0, 0.0),
                        (acc, body) -> acc.accumForce(body),
                        (acc1, acc2) -> acc1.merge(acc2),
                        acc -> new Body(b, acc.x, acc.y));
        return bodies.stream()
                .map(body -> b.calculateForce(body))
                .collect(forceCollector);
    }

    public static ArrayList<Body> moveBodies(List<Body> bodies){
        return bodies.stream()
                .map(NBody::move)
                .collect(Collectors.toCollection(ArrayList::new));
    }
    public static Body move(Body b) {
        Body ret = new Body(b);
        double dvX, dvY, dpX, dpY;
        dvX = ret.fx / ret.mass * NBody.DT;
        dvY = ret.fy / ret.mass * NBody.DT;
        dpX = (ret.vx + dvX/2.0) * NBody.DT;
        dpY = (ret.vy + dvY/2.0) * NBody.DT;

        ret.vx += dvX;
        ret.vy += dvY;
        ret.setLocation(b.getX() + dpX, b.getY() + dpY);
        ret.fx = 0.0; //reset force vector?
        ret.fy = 0.0;
        //System.out.print(String.format("(%.2f,%.2f) ",ret.getX(),ret.getY()));
        return ret;
    }

    public static ArrayList<Body> checkCollisions(List<Body> bodies){
        return bodies.stream()
                .map(b -> NBody.bodies.stream() //calculate against all
                        .filter(b2 -> hasCollided(b,b2))
                        .reduce(b, (b1, b2) -> b1.checkCollision(b2)))
                .map(NBody::checkWall)
                .collect(Collectors.toCollection(ArrayList::new));

    }
    private static Body checkWall(Body b){
        Body ret = new Body(b);
        if(b.hitWall()) {
            ret.vx = -ret.vx;
            ret.vy = -ret.vy;
        }
        return ret;
    }

    public static boolean hasCollided(Body body1, Body body2){
        return body1.distance(body2) < (body1.radius + body2.radius);
    }

    public static void barrier(int pid) throws InterruptedException {
        for(int i = 0; i < numStages; i++) {
            int partner = (pid + ipow(2, i)) % numWorkers;
            stages[pid][i].release();
            stages[partner][i].acquire();
        }
    }

    private static int ipow(int base, int exp) {
        int result = 1;
        while (exp != 0) {
            if((exp&1) != 0) result *= base;
            exp >>= 1;
            base *= base;
        }
        return result;
    }

}
class Accum {
    double x, y;

    public Accum(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public Accum accumForce(Body b) {
        this.x += b.fx;
        this.y += b.fy;
        return new Accum(this.x, this.y);
    }

    public Accum merge(Accum acc){
        this.x += acc.x;
        this.y += acc.y;
        return new Accum(this.x, this.y);
    }
}
class Body extends Point2D.Double {
    double vx, vy; //velocity vector
    double mass, radius; //mass in somehrting
    double fx, fy; //acting force
    String label;

    public Body(double x, double y, double radius, double mass) {
        this.setLocation(x, y);
        this.radius = radius;
        this.mass = mass;
        this.makeLabel();

    }

    private void makeLabel() {
        label = String.format("(%.2f, %.2f)\nr = %.0f\nmass = %.0f", this.getX(), this.getY(), radius, mass);
    }

    public Body(Body b) {
        this.setLocation(b.getX(), b.getY());
        this.radius = b.radius;
        this.mass = b.mass;
        this.fx = b.fx;
        this.fy = b.fy;
        this.vx = b.vx;
        this.vy = b.vy;
        this.makeLabel();
    }

    public Body(Body b, double fx, double fy) { // for accumulating force
        this.setLocation(b.getX(), b.getY());
        this.radius = b.radius;
        this.mass = b.mass;
        this.fx = fx;
        this.fy = fy;
        this.vx = b.vx;
        this.vy = b.vy;
        this.makeLabel();
    }

    public boolean hitWall() {
        return getX() - radius < 0 ||
                getX() + radius > NBody.framesize ||
                getY() - radius < 0 ||
                getY() + radius > NBody.framesize;
    }

    public Body calculateForce(Body b) {
        double distance, magnitude, directionX, directionY, forceX, forceY;
        distance = this.distance(b);
        if (distance == 0) return this;

        magnitude = (NBody.G * this.mass * b.mass) / Math.pow(distance, 2);
//        System.out.println(String.format("Magnitude is %f",magnitude));
        directionX = b.getX() - this.getX();
        directionY = b.getY() - this.getY();
        forceX = magnitude * directionX / distance;
        forceY = magnitude * directionY / distance;

        return new Body(this, forceX, forceY);
    }

    public Body checkCollision(Body b) {
        //no one is going to have a radius smaller than this, right?
        if (this.equals(b) || (Math.abs(this.getX() - b.getX()) <= .1 && Math.abs(this.getY() - b.getY()) <= .1))
            return new Body(this);
        double newVx, newVy;
        synchronized(NBody.collisionSync) {
            NBody.collisions++;
        }
        newVx = (this.vx * (this.mass - b.mass) +
                (2.0 * b.mass * b.vx)) / (this.mass + b.mass);
        newVy = (this.vy * (this.mass - b.mass) +
                (2.0 * b.mass * b.vy)) / (this.mass + b.mass);
//        System.out.println(String.format("Collision between bodies at (%.2f, %.2f) and (%.2f, %.2f)",
 //               this.getX(), this.getY(), b.getX(), b.getY()));
        Body ret = new Body(this);
        ret.vx = newVx;
        ret.vy = newVy;
        return ret;
    }
}

class NBodyWorker extends Thread {
    //something like id*n -> id*n + n
    //for assigned bodies
    int id, start, end, steps; //right exclusive
    public List<Body> bodies;
    public NBodyWorker(int id, int start, int end, int steps) {
        this.id = id;
        this.start = start;
        this.end = end;
        this.steps = steps;
    }
    @Override
    public void run() {
        try {
            NBody.barrier(id);
            for (int i = 0; i < steps; i++) {
                bodies = NBody.bodies.subList(start, end);
                bodies = NBody.moveBodies(NBody.checkCollisions(NBody.calculateForces(bodies)));
                NBody.barrier(id);
                if(id == 0) {
                    NBody.forMain.release();
                    NBody.forWorker.acquire();
                }
                NBody.barrier(id);
            }
        }
        catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
