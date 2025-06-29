#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <sstream>
#include <map>

using namespace std;
namespace fs = std::filesystem;

struct GeometryConfig {
    double X_end;
    double Y_end;
    int X_direction;
    int Y_direction;
};

struct BoundaryConfig {
    double top_velocity;
    double reynolds_number;
};

struct InitialConfig {
    double initial_top_velocity;
};

struct SolverConfig {
    std::string poisson_scheme;
    std::string vorticity_scheme;
    double delta_t;
};

// Utility Functions
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::string removeQuotes(const std::string& str) {
    std::string trimmed = trim(str);
    if ((trimmed.front() == '"' || trimmed.front() == '\'') && trimmed.back() == trimmed.front())
        return trimmed.substr(1, trimmed.size() - 2);
    return trimmed;
}

// Functions to read data from files
bool readGeometryConfig(const std::string& filename, GeometryConfig& config) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line, current_section;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line.substr(0, 2) == "//") continue;

        if (line.find("Domain_info") != std::string::npos) {
            current_section = "Domain_info"; continue;
        } else if (line.find("Grid_points") != std::string::npos) {
            current_section = "Grid_points"; continue;
        }

        if (line == "{" || line == "}") continue;

        size_t eq_pos = line.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = trim(line.substr(0, eq_pos));
            std::string value = trim(line.substr(eq_pos + 1));
            if (!value.empty() && value.back() == ';') value.pop_back();

            if (current_section == "Domain_info") {
                if (key == "X_end") config.X_end = std::stod(value);
                else if (key == "Y_end") config.Y_end = std::stod(value);
            } else if (current_section == "Grid_points") {
                if (key == "X_direction") config.X_direction = std::stoi(value);
                else if (key == "Y_direction") config.Y_direction = std::stoi(value);
            }
        }
    }

    file.close();
    return true;
}

bool readBoundaryConfig(const std::string& filename, BoundaryConfig& config) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line, current_section;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line.substr(0, 2) == "//") continue;

        if (line.find("U_velocity") != std::string::npos) {
            current_section = "U_velocity"; continue;
        } else if (line.find("Reynolds_Number") != std::string::npos) {
            current_section = "Reynolds_Number"; continue;
        }

        if (line == "{" || line == "}") continue;

        size_t eq_pos = line.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = trim(line.substr(0, eq_pos));
            std::string value = trim(line.substr(eq_pos + 1));
            if (!value.empty() && value.back() == ';') value.pop_back();

            if (current_section == "U_velocity" && key == "top")
                config.top_velocity = std::stod(value);
            else if (current_section == "Reynolds_Number" && key == "Re")
                config.reynolds_number = std::stod(value);
        }
    }

    file.close();
    return true;
}

bool readInitialConfig(const std::string& filename, InitialConfig& config) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    bool in_u_velocity = false;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line.substr(0, 2) == "//") continue;

        if (line.find("U_velocity") != std::string::npos) {
            in_u_velocity = true; continue;
        }

        if (line == "{") continue;
        if (line == "}") { in_u_velocity = false; continue; }

        if (in_u_velocity) {
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = trim(line.substr(0, eq_pos));
                std::string value = trim(line.substr(eq_pos + 1));
                if (!value.empty() && value.back() == ';') value.pop_back();
                if (key == "top") config.initial_top_velocity = std::stod(value);
            }
        }
    }

    file.close();
    return true;
}

bool readSolverConfig(const std::string& filename, SolverConfig& config) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line, current_section;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line.substr(0, 2) == "//") continue;

        if (line.find("Poisson") != std::string::npos) {
            current_section = "Poisson"; continue;
        } else if (line.find("Vorticity") != std::string::npos) {
            current_section = "Vorticity"; continue;
        } else if (line.find("Timestep") != std::string::npos) {
            current_section = "Timestep"; continue;
        }

        if (line == "{" || line == "}") continue;

        size_t eq_pos = line.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = trim(line.substr(0, eq_pos));
            std::string value = trim(line.substr(eq_pos + 1));
            if (!value.empty() && value.back() == ';') value.pop_back();
            value = removeQuotes(value);

            if (current_section == "Poisson" && key == "Scheme")
                config.poisson_scheme = value;
            else if (current_section == "Vorticity" && key == "Scheme")
                config.vorticity_scheme = value;
            else if (current_section == "Timestep" && key == "delta_t")
                config.delta_t = std::stod(value);
        }
    }

    file.close();
    return true;
}

// Unified Reader
bool readAllConfigs(GeometryConfig& geom, BoundaryConfig& boundary, InitialConfig& initial, SolverConfig& solver) {
    return readGeometryConfig("../geometry/geom.txt", geom)
        && readBoundaryConfig("../boundaryConditions/boundaryCond.txt", boundary)
        && readInitialConfig("../initialConditions/initialCond.txt", initial)
        && readSolverConfig("../solverSettings/solverProperties.txt", solver);
}

// The following function initializes the values of the velocities, stream function and vorticity to the domain
void initialize_domain(int gridpts_x, int gridpts_y, double delta_x, double delta_y, double initial_top_lid_velocity, vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& shi, vector<vector<double>>& zeta)
{
    // u and v initialization to domain
    for(int j=0;j<gridpts_y;j++)
    {
        if(j==0)
        {
            for(int i=0;i<gridpts_x;i++)
            {
                u[j][i] = initial_top_lid_velocity;
                v[j][i] = 0;
            }
        }
        else
        {
            for(int i=0;i<gridpts_x;i++)
            {
                u[j][i] = 0;
                v[j][i] = 0;
            }
        }
    }

    // zeta (vorticity) initialization
    for(int j=1;j<gridpts_y-1;j++)
    {
        for(int i=1;i<gridpts_x-1;i++)
        {
            zeta[j][i] = (v[j][i+1] - v[j][i-1])/(2*delta_x) - (u[j+1][i] - u[j-1][i])/(2*delta_y);
        }
    }

    for(int i=1;i<gridpts_x-1;i++)
    {
        zeta[0][i] = (v[0][i+1] - v[0][i-1])/(2*delta_x) - (u[1][i] - u[0][i])/(delta_y);
        zeta[gridpts_y-1][i] = (v[gridpts_y-1][i+1] - v[gridpts_y-1][i-1])/(2*delta_x) - (u[gridpts_y-1][i] - u[gridpts_y-2][i])/(delta_y);
    }

    for(int j=1;j<gridpts_y-1;j++)
    {
        zeta[j][0] = (v[j][1] - v[j][0])/(delta_x) - (u[j+1][0] - u[j-1][0])/(2*delta_y);
        zeta[j][gridpts_x-1] = (v[j][gridpts_x-1] - v[j][gridpts_x-2])/(delta_x) - (u[j+1][gridpts_x-1] - u[j-1][gridpts_x-1])/(2*delta_y);
    }

    zeta[0][0] = (v[0][1] - v[0][0])/delta_x - (u[1][0] - u[0][0])/delta_y;
    zeta[0][gridpts_x-1] = (v[0][gridpts_x-1] - v[0][gridpts_x-2])/delta_x - (u[1][gridpts_x-1] - u[0][gridpts_x-1])/delta_y;
    zeta[gridpts_y-1][0] = (v[gridpts_y-1][1] - v[gridpts_y-1][0])/delta_x - (u[1][0] - u[0][0])/delta_y;
    zeta[gridpts_y-1][gridpts_x-1] = (v[gridpts_y-1][gridpts_x-1] - v[gridpts_y-1][gridpts_x-2])/delta_x - (u[gridpts_y-1][gridpts_x-1] - u[gridpts_y-2][gridpts_x-1])/delta_y;

    // Shi (Stream Function) initialization
    for(int j=0;j<gridpts_y;j++)
    {
        for(int i=0;i<gridpts_x;i++)
        {
            shi[j][i] = 0;
        }
    }
}

// The following function marches the vorticity values for the internal grid points to the next timestep
void FTCS(int gridpts_x, int gridpts_y, double delta_x, double delta_y, double delta_t, double Re, vector<vector<double>>& zeta, vector<vector<double>>& copy_zeta, vector<vector<double>>& u, vector<vector<double>>& v)
{
    for(int j=0;j<gridpts_y;j++)
    {
        for(int i=0;i<gridpts_x;i++)
        {
            copy_zeta[j][i] = zeta[j][i];
        }
    }

    // FTCS Scheme implemented for the vorticity equation
    for(int j=1;j<gridpts_y-1;j++)
    {
        for(int i=1;i<gridpts_x-1;i++)
        {
            double term1 = copy_zeta[j][i] - u[j][i]*(delta_t/(2*delta_x))*(copy_zeta[j][i+1] - copy_zeta[j][i-1]) - v[j][i]*(delta_t/(2*delta_y))*(copy_zeta[j+1][i] - copy_zeta[j-1][i]);
            double term2 = (delta_t/Re)*((copy_zeta[j][i-1] - 2*copy_zeta[j][i] + copy_zeta[j][i+1])/pow(delta_x,2) + (copy_zeta[j-1][i] - 2*copy_zeta[j][i] + copy_zeta[j+1][i])/pow(delta_y,2));
            zeta[j][i] = term1 + term2;
        }
    }
}

// This function evaluates the residual value in the stream function values
double shi_convergence_check_val(int gridpts_x, int gridpts_y, vector<vector<double>>& shi, vector<vector<double>>& copy_shi)
{
    double shi_sum = 0;
    for(int j=0;j<gridpts_y;j++)
    {
        for(int i=0;i<gridpts_x;i++)
        {
            shi_sum = shi_sum + pow(shi[j][i] - copy_shi[j][i],2);
        }
    }
    return sqrt(shi_sum);
}

// The following function solves the Poisson equation using the SOR (Successive Over Relaxation) Scheme into the Gauss-Siedel Update value
void GS_with_SOR(int gridpts_x, int gridpts_y, double delta_x, double delta_y, double shi_tolerance, double omega, vector<vector<double>>& shi, vector<vector<double>>& copy_shi, vector<vector<double>>& zeta)
{
    double beta = delta_x/delta_y;

    do
    {
        for(int j=0;j<gridpts_y;j++)
        {
            for(int i=0;i<gridpts_x;i++)
            {
                copy_shi[j][i] = shi[j][i];
            }
        }

        for(int j=1;j<gridpts_y-1;j++)
        {
            for(int i=1;i<gridpts_x-1;i++)
            {
                double shi_gs_update = (copy_shi[j][i+1] + shi[j][i-1] + copy_shi[j+1][i] + shi[j-1][i] + zeta[j][i]*pow(delta_x,2))/(2*(1+beta*beta));
                shi[j][i] = (1 - omega)*copy_shi[j][i] + omega*shi_gs_update;
            }
        }
    }
    while(shi_convergence_check_val(gridpts_x,gridpts_y,shi,copy_shi) >= shi_tolerance);
}

// The following function evaluates the velocity component values from the stream function values
void update_velocities(int gridpts_x, int gridpts_y, double delta_x, double delta_y, double boundary_top_lid_velocity, vector<vector<double>>& shi, vector<vector<double>>& u, vector<vector<double>>& v)
{
    for(int j=1;j<gridpts_y-1;j++)
    {
        for(int i=1;i<gridpts_x-1;i++)
        {
            u[j][i] = (shi[j+1][i] - shi[j-1][i])/(2*delta_y);
            v[j][i] = -(shi[j][i+1] - shi[j][i-1])/(2*delta_x);
        }
    }

    for(int i=0;i<gridpts_x;i++)
    {
        u[0][i] = boundary_top_lid_velocity;
    }
}

// The following function updates the boundary vorticity values of the domain
void update_zeta_boundaries(int gridpts_x, int gridpts_y, double delta_x, double delta_y, vector<vector<double>>& zeta, vector<vector<double>>& shi, vector<vector<double>>& u)
{
    for(int i=1;i<gridpts_x-1;i++)
    {
        zeta[0][i] = (2/pow(delta_y,2))*(u[0][i]*delta_y - shi[1][i]);
        zeta[gridpts_y-1][i] = (2/pow(delta_y,2))*(-shi[gridpts_y-2][i]);
    }

    for(int j=1;j<gridpts_y-1;j++)
    {
        zeta[j][0] = (2/pow(delta_x,2))*(-shi[j][1]);
        zeta[j][gridpts_x-1] = (2/pow(delta_x,2))*(-shi[j][gridpts_x-2]);
    }
}

// This function evaluates the residual in the vorticity value
double zeta_convergence_check_val(int gridpts_x, int gridpts_y, vector<vector<double>>& zeta, vector<vector<double>>& copy_zeta)
{
    double zeta_sum = 0;
    for(int j=0;j<gridpts_y;j++)
    {
        for(int i=0;i<gridpts_x;i++)
        {
            zeta_sum = zeta_sum + pow(zeta[j][i] - copy_zeta[j][i],2);
        }
    }
    return sqrt(zeta_sum);
}

/* The following function inverts the sign of the required physical quantity due to y-axis being taken in an inverted configuration
   The values of v (vertical velocity), stream function and vorticity need to be inverted because of that. */
void invert_sign_due_to_axes_convection(int gridpts_x, int gridpts_y, vector<vector<double>>& variable)
{
    for(int j=0;j<gridpts_y;j++)
    {
        for(int i=0;i<gridpts_x;i++)
        {
            variable[j][i] = -variable[j][i];
        }
    }
}

// Writes domain information in CSV file
void write_domain(const string& filename, int gridpts_x, int gridpts_y, double delta_x, double delta_y)
{
    ofstream file(filename);
    double x_domain_counter = 0, y_domain_counter = 0;

    for(int i=0;i<gridpts_x;i++)
    {
        file << x_domain_counter;
        if(i<gridpts_x-1)
        {
            file << ",";
        }
        x_domain_counter = x_domain_counter + delta_x;
    }
    file << "\n";

    for(int j=0;j<gridpts_y;j++)
    {
        file << y_domain_counter;
        if(j<gridpts_y-1)
        {
            file << ",";
        }
        y_domain_counter = y_domain_counter + delta_y;
    }
    file << endl;
    file.close();
}

// Writes output values in CSV file
void write_output_in_csv(const string& filename, int gridpts_x, int gridpts_y, vector<vector<double>>& variable)
{
    ofstream file(filename);
    for(int j=0;j<gridpts_y;j++)
    {
        for(int i=0;i<gridpts_x;i++)
        {
            file << variable[gridpts_y-1-j][i];
            if(i<gridpts_x-1)
            {
                file << ",";
            }
        }
        if (j<gridpts_y-1)
        {
            file << "\n";
        }
    }
    file << endl;
    file.close();
}

int main()
{
    // Variable names are self explanatory
    GeometryConfig geom;
    BoundaryConfig boundary;
    InitialConfig initial;
    SolverConfig solver;

    if (!readAllConfigs(geom, boundary, initial, solver)) {
        cerr << "Error: Failed to load configuration files.\n";
        return 1;
    }

    cout << "Configuration loaded successfully:\n";
    cout << "Domain: " << geom.X_end << " x " << geom.Y_end << "\n";
    cout << "Grid: " << geom.X_direction << " x " << geom.Y_direction << "\n";
    cout << "Initial top velocity: " << initial.initial_top_velocity << "\n";
    cout << "Boundary top velocity: " << boundary.top_velocity << "\n";
    cout << "Reynolds number: " << boundary.reynolds_number << "\n";
    cout << "Poisson scheme: " << solver.poisson_scheme << "\n";
    cout << "Vorticity scheme: " << solver.vorticity_scheme << "\n";
    cout << "Time step: " << solver.delta_t << "\n\n";

    double domain_length_x = geom.X_end;
    double domain_length_y = geom.Y_end;
    int gridpts_x = geom.X_direction;
    int gridpts_y = geom.Y_direction;
    double initial_top_lid_velocity = initial.initial_top_velocity;
    double boundary_top_lid_velocity = boundary.top_velocity;
    double Re = boundary.reynolds_number;
    double delta_t = solver.delta_t;

    double delta_x = domain_length_x / (gridpts_x - 1);
    double delta_y = domain_length_y / (gridpts_y - 1);

    ostringstream ss_gptsx, ss_gptsy, ss_Re;
    ss_gptsx << fixed << setprecision(0) << gridpts_x;
    ss_gptsy << fixed << setprecision(0) << gridpts_y;
    ss_Re    << fixed << setprecision(0) << Re;

    string gridpts_x_string = ss_gptsx.str();
    string gridpts_y_string = ss_gptsy.str();
    string Re_string = ss_Re.str();

    string poisson_solver_name = solver.poisson_scheme;
    string vorticity_scheme_name = solver.vorticity_scheme;

    double shi_tolerance = pow(10,-3);
    double zeta_tolerance = pow(10,-3);
    double omega = 2.0/(1 + sqrt(1 - (cos(M_PI/(gridpts_x-1)) + cos(M_PI/(gridpts_y-1)))/2));        // SOR (Successive Over Relaxation) factor

    vector<vector<double>> u(gridpts_y, vector<double>(gridpts_x));
    vector<vector<double>> v(gridpts_y, vector<double>(gridpts_x));
    vector<vector<double>> shi(gridpts_y, vector<double>(gridpts_x));
    vector<vector<double>> zeta(gridpts_y, vector<double>(gridpts_x));
    vector<vector<double>> copy_shi(gridpts_y, vector<double>(gridpts_x));
    vector<vector<double>> copy_zeta(gridpts_y, vector<double>(gridpts_x));

    double time_counter = 0;
    int output_counter = 0;
    double interval_after_which_store_result = 0.5;
    int output_interval = (int)(interval_after_which_store_result/delta_t);
    int global_counter = 0;

    initialize_domain(gridpts_x,gridpts_y,delta_x,delta_y,initial_top_lid_velocity,u,v,shi,zeta);

    do
    {
        if(global_counter == 0)
        {
            ostringstream ss_first_timestep;
            ss_first_timestep << fixed << setprecision(2) << time_counter;
            string time_string = ss_first_timestep.str();
            string output_dir = string("../postProcessing/") + "Test_Case_Re_" + Re_string + "_Grid_" + gridpts_x_string + "x" + gridpts_y_string + "/" + time_string;
            fs::create_directories(output_dir);
            write_domain(output_dir + "/Domain.csv", gridpts_x, gridpts_y, delta_x, delta_y);
            write_output_in_csv(output_dir + "/StreamFunction.csv", gridpts_x, gridpts_y, shi);
            write_output_in_csv(output_dir + "/Vorticity.csv", gridpts_x, gridpts_y, zeta);
            write_output_in_csv(output_dir + "/U.csv", gridpts_x, gridpts_y, u);
            write_output_in_csv(output_dir + "/V.csv", gridpts_x, gridpts_y, v);
        }

        if(vorticity_scheme_name == "FTCS")
        {
            FTCS(gridpts_x,gridpts_y,delta_x,delta_y,delta_t,Re,zeta,copy_zeta,u,v);
        }

        if(poisson_solver_name == "GS_with_SOR")
        {
            GS_with_SOR(gridpts_x,gridpts_y,delta_x,delta_y,shi_tolerance,omega,shi,copy_shi,zeta);
        }

        update_velocities(gridpts_x,gridpts_y,delta_x,delta_y,boundary_top_lid_velocity,shi,u,v);
        update_zeta_boundaries(gridpts_x,gridpts_y,delta_x,delta_y,zeta,shi,u);

        time_counter = time_counter + delta_t;
        output_counter++;
        global_counter++;

        if(output_counter >= output_interval)
        {
            output_counter = 0;
            ostringstream ss_time;
            ss_time << fixed << setprecision(2) << time_counter;
            string time_string = ss_time.str();

            string output_dir = string("../postProcessing/") + "Test_Case_Re_" + Re_string + "_Grid_" + gridpts_x_string + "x" + gridpts_y_string + "/" + time_string;
            fs::create_directories(output_dir);

            invert_sign_due_to_axes_convection(gridpts_x,gridpts_y,v);
            invert_sign_due_to_axes_convection(gridpts_x,gridpts_y,shi);
            invert_sign_due_to_axes_convection(gridpts_x,gridpts_y,zeta);

            //write_domain(output_dir + "/Domain.csv", gridpts_x, gridpts_y, delta_x, delta_y);
            write_output_in_csv(output_dir + "/StreamFunction.csv", gridpts_x, gridpts_y, shi);
            write_output_in_csv(output_dir + "/Vorticity.csv", gridpts_x, gridpts_y, zeta);
            write_output_in_csv(output_dir + "/U.csv", gridpts_x, gridpts_y, u);
            write_output_in_csv(output_dir + "/V.csv", gridpts_x, gridpts_y, v);

            invert_sign_due_to_axes_convection(gridpts_x,gridpts_y,v);
            invert_sign_due_to_axes_convection(gridpts_x,gridpts_y,shi);
            invert_sign_due_to_axes_convection(gridpts_x,gridpts_y,zeta);

            cout << "Simulation complete till: " << time_counter << " sec\n" << "Residue = " << zeta_convergence_check_val(gridpts_x,gridpts_y,zeta,copy_zeta) << "\n";
        }
    }
    while(zeta_convergence_check_val(gridpts_x,gridpts_y,zeta,copy_zeta) >= zeta_tolerance);

    cout << "Successful convergence achieved! \nTime taken for steady state solution is: " << time_counter << " sec\n\n";

    string case_folder = "Test_Case_Re_" + Re_string + "_Grid_" + gridpts_x_string + "x" + gridpts_y_string;

    // MATLAB call commmand
    string matlab_command = "cd ../executables && ";
    matlab_command += "set MATLAB_CASE_NAME=" + case_folder + " && ";
    matlab_command += "matlab -batch run_make_variable_video";

    // Execution of MATLAB call command
    int ret = system(matlab_command.c_str());
    if (ret != 0)
    {
        cerr << "MATLAB script execution failed\n";
    }
    else
    {
        cout << "MATLAB script executed successfull\n";
    }
}
