{
  description = "Astrocalibrator - A lightweight calibration tool for astrophotography workflows";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        
        pythonEnv = pkgs.python311.withPackages (ps: with ps; [
          # Core dependencies
          pillow
          astropy
          astroquery
          scipy
          numpy
          
          # Tkinter is required for the GUI
          tkinter
          
          # Development tools
          ipython
          black
          pylint
        ]);
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pythonEnv
            pkgs.git
          ];
          
          shellHook = ''
            echo "Welcome to Astrocalibrator development environment!"
            echo "Python version: $(python --version)"
            echo ""
            echo "To run the application:"
            echo "  python main.py"
            echo ""
            echo "Note: You'll need to install ASTAP separately for plate solving:"
            echo "  Download from: https://www.hnsky.org/astap.htm"
            echo "  Configure the ASTAP executable location in the app settings"
          '';
          
          # Set environment variables if needed
          PYTHONPATH = "./";
        };
      }
    );
}
