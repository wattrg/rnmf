extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn;

#[proc_macro_derive(UserConfig)]
pub fn user_config_derive(input: TokenStream) -> TokenStream {
    // Construct a representation of Rust code as a syntax tree
    // that we can manipulate
    let ast = syn::parse(input).unwrap();

    // Build the trait implementation
    impl_user_config_macro(&ast)
}

fn impl_user_config_macro(ast: &syn::DeriveInput)->TokenStream{
    let name = &ast.ident;
}