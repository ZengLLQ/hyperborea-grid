use clap::Parser;
use hyperborea_grid::cli::Cli;

fn main() {
    let cli = Cli::parse();
    cli.execute();
}
